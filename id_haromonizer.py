#!/usr/bin/env python3
"""ID Haromonizer.

Workflow:
  1. Reads CSV file(s) from a local `csv/` folder.
  2. Prompts the user to select the column that holds the patient id
     (patient id / mr_number - column name may be arbitrary).
  3. Copies the CSV to the output folder with the patient id column
     stripped of alphabets and with all special characters / spaces
     converted to underscores.
  4. Reads all DICOM files found inside zip files (nested zips supported)
     from a local `dicom/` folder and extracts the patient ids.
  5. Matches DICOM patient ids to the cleaned CSV patient ids and copies
     the matching DICOM files to the output folder organized by
     Study > Series > Instance, file name being
     <SOP Instance UID>_<Modality>.dcm (e.g. ...1229678._CT.dcm).
     (0010,0020) PatientID inside each written file is set to the cleaned id the
     files were matched on (files already carrying it are copied byte for byte);
     no other tag is touched.
     One file per SOP Instance UID: if the same instance is present twice
     (same UID in several members/archives - on this corpus only
     (0008,0013) InstanceCreationTime differs between such copies) the copy
     that comes last in archive order wins, i.e. duplicates overwrite.
  6. Unmatched DICOM files are kept in a separate unmatched folder.
  7. Writes result CSVs:
       - patients with no matching DICOM ("not found")
       - patients with matching DICOM (patient, studies, series, deduped
         modalities, instance count, folder paths)
       - unmatched DICOM patients (patient, instance count, modalities and the
         source archive each came from) as results/unmatched.csv

The original input CSV files are never modified.

Memory / IO model
------------------
The script never holds DICOM payloads in RAM. It runs in two sequential passes:

  pass 1  opens each archive for *random access* and parses DICOM tags only
          (pixel data is skipped), keeping a small metadata index of
          {patient id, UIDs, archive, member name} per instance;
  pass 2  reopens each archive once and streams the needed members straight to
          disk with `shutil.copyfileobj` (1 MiB chunks), except for files whose
          (0010,0020) PatientID differs from the id it was matched on: those are
          re-read in full and resaved with the corrected id (one file in RAM).

Peak RAM is therefore a few MB plus the metadata index, independent of the
corpus size, and progress is printed continuously so the run is never silent.
"""

import argparse
import contextlib
import csv
import hashlib
import logging
import os
import re
import shutil
import tempfile
import time
import warnings
import zipfile
from pathlib import Path

import pydicom

DEFAULT_CSV_DIR = "csv"
DEFAULT_DICOM_DIR = "dicom"
DEFAULT_OUTPUT_DIR = "output"

COPY_CHUNK = 1 << 20               # streaming buffer for every byte copy (1 MiB)
SPOOL_LIMIT = 64 * 1024 * 1024     # nested zips bigger than this spill to disk

# pydicom is very chatty about non-conformant files; we handle those quietly and
# report a single skipped-file count per archive instead.
logging.getLogger("pydicom").setLevel(logging.ERROR)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def normalize_id(value) -> str:
    """Strip alphabets and convert all special characters and spaces to '_'.

    Examples:
        "MRN12345"      -> "12345"
        "MRN-12345"     -> "12345"
        "MRN 123/45"    -> "123_45"
        "PatientABC123" -> "123"
    """
    s = str(value).strip().upper()
    # remove alphabets (A-Z after upper-casing)
    s = re.sub(r"[A-Z]+", "", s)
    # convert any remaining special character / space to a single underscore
    s = re.sub(r"[^0-9]+", "_", s)
    # trim leading/trailing underscores
    s = s.strip("_")
    return s


def _elem_value(ds, keyword, default=""):
    """Safely read a string value from a pydicom dataset by keyword."""
    try:
        elem = ds.data_element(keyword)
    except Exception:
        return default
    if elem is None or elem.value is None:
        return default
    return str(elem.value).strip()


def _safe_component(value: str) -> str:
    """Make a UID safe to use as a path component (force-mode files can be odd)."""
    value = re.sub(r"[^0-9.]+", "_", str(value)).strip("._")
    return value or "UNKNOWN"


def _modality_token(value: str) -> str:
    """Modality as a safe filename token: 'CT', 'RTSTRUCT', 'MR'... (max 16 chars)."""
    token = re.sub(r"[^0-9A-Za-z]+", "", str(value)).upper()[:16]
    return token or "UNKNOWN"


# --------------------------------------------------------------------------
# Pass 1 - metadata only
# --------------------------------------------------------------------------
@contextlib.contextmanager
def _spool_member(zf, info):
    """Stream a zip member to a spooled file (RAM up to SPOOL_LIMIT, then disk).

    Used for *nested* zips only: it gives random access to the inner archive
    without ever holding the whole thing in memory.
    """
    tmp = tempfile.SpooledTemporaryFile(max_size=SPOOL_LIMIT)
    try:
        with zf.open(info) as src:
            shutil.copyfileobj(src, tmp, COPY_CHUNK)
        tmp.seek(0)
        yield tmp
    finally:
        tmp.close()


def _iter_entries(zf, label, nested, processed):
    """Depth-first walk of an open archive; yields (label, nested, zf, info).

    `nested` is the chain of member names leading from the top-level archive to
    the archive that holds the entry (empty for plain archives) so pass 2 can
    find the bytes again without anything being kept open or in memory.
    """
    for info in zf.infolist():
        if info.is_dir():
            continue
        name = info.filename
        if name.lower().endswith(".zip"):
            # Dedupe nested archives by content signature taken straight from the
            # central directory (size + CRC) - no decompression, no hashing.
            sig = (info.file_size, info.CRC)
            if sig in processed:
                continue
            processed.add(sig)
            with _spool_member(zf, info) as fh:
                try:
                    inner = zipfile.ZipFile(fh)
                except (zipfile.BadZipFile, OSError):
                    continue
                try:
                    yield from _iter_entries(
                        inner,
                        f"{label}!{name}",
                        nested + (name,),
                        processed,
                    )
                finally:
                    inner.close()
        else:
            yield f"{label}!{name}", nested, zf, info


def _read_header(zf, info):
    """Read *tags only* from a zip member - pixel data is never loaded."""
    for kwargs in ({}, {"force": True}):
        try:
            with zf.open(info) as raw:
                return pydicom.dcmread(raw, stop_before_pixels=True, **kwargs)
        except Exception:
            continue
    return None


def _record(ds, source, zip_path, nested, member, size):
    pid = _elem_value(ds, "PatientID") or _elem_value(ds, "PatientName")
    sop_uid = _safe_component(_elem_value(ds, "SOPInstanceUID"))
    modality = _modality_token(_elem_value(ds, "Modality"))
    return {
        "patient_id": pid,
        "cleaned_id": normalize_id(pid),
        "study_uid": _safe_component(_elem_value(ds, "StudyInstanceUID")),
        "series_uid": _safe_component(_elem_value(ds, "SeriesInstanceUID")),
        "sop_uid": sop_uid,
        "modality": modality,
        # <SOP Instance UID>_<Modality>.dcm
        "instance_name": f"{sop_uid}_{modality}.dcm",
        # Where the bytes live: reopened lazily in pass 2, never kept in RAM.
        "zip_path": Path(zip_path),
        "nested": tuple(nested),
        "member": member,
        "size": size,
        "source": source,
    }


def _zip_signature(path: Path):
    """Cheap identity for an archive: central directory only, no data read."""
    with zipfile.ZipFile(path, "r") as zf:
        infos = zf.infolist()
        return (
            path.stat().st_size,
            len(infos),
            sum(i.file_size for i in infos),
            sum(i.CRC for i in infos),
        )


def _hash_file(path: Path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(COPY_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def build_dicom_index(dicom_dir: Path, progress=True):
    """Pass 1: index DICOM *tags* by normalized patient id.

    Only tag values and zip coordinates are retained, so memory stays flat no
    matter how large the corpus is (a 12 GB corpus costs a few MB of RAM).
    """
    index = {}
    dicom_dir = Path(dicom_dir)
    if not dicom_dir.exists():
        print(f"DICOM folder '{dicom_dir}' does not exist.")
        return index

    zip_paths = sorted(p for p in dicom_dir.rglob("*.zip") if p.is_file())
    print(f"\nPass 1/2: reading DICOM headers from {len(zip_paths)} archive(s)")

    nested_seen = set()
    sigs = {}
    total_recs = 0
    started = time.time()

    for n, zf_path in enumerate(zip_paths, 1):
        try:
            sig = _zip_signature(zf_path)
        except (zipfile.BadZipFile, OSError) as exc:
            print(f"  [{n}/{len(zip_paths)}] {zf_path.name}: not a readable zip ({exc}) - skipped")
            continue

        if sig in sigs:
            # Identical signature: confirm with a real hash before discarding.
            if _hash_file(sigs[sig]) == _hash_file(zf_path):
                print(f"  [{n}/{len(zip_paths)}] {zf_path.name}: duplicate of {sigs[sig].name} - skipped")
                continue
        sigs[sig] = zf_path

        count = skipped = 0
        local_pids = set()
        try:
            with zipfile.ZipFile(zf_path, "r") as zf:
                for source, nested, live_zf, info in _iter_entries(zf, zf_path.name, (), nested_seen):
                    ds = _read_header(live_zf, info)
                    if ds is None:
                        skipped += 1
                        continue
                    rec = _record(ds, source, zf_path, nested, info.filename, info.file_size)
                    if not rec["cleaned_id"]:
                        skipped += 1
                        continue
                    index.setdefault(rec["cleaned_id"], []).append(rec)
                    local_pids.add(rec["cleaned_id"])
                    count += 1
        except (zipfile.BadZipFile, OSError) as exc:
            print(f"  [{n}/{len(zip_paths)}] {zf_path.name}: read error ({exc}) - partially indexed")

        total_recs += count
        if progress:
            print(f"  [{n}/{len(zip_paths)}] {zf_path.name}: {count} DICOM files"
                  f" for {len(local_pids)} patients"
                  f"{f', {skipped} unreadable' if skipped else ''}"
                  f" ({total_recs} total, {time.time() - started:.0f}s)")

    for pid in index:
        index[pid].sort(key=lambda r: (r["study_uid"], r["series_uid"], r["sop_uid"]))
    print(f"  Indexed {total_recs} DICOM files for {len(index)} patients"
          f" in {time.time() - started:.0f}s")
    return index


# --------------------------------------------------------------------------
# Pass 2 - stream matched / unmatched files to disk
# --------------------------------------------------------------------------
class MemberReader:
    """Streams members out of the source archives, one archive open at a time.

    Records are fed in archive order, so each archive (and each nested archive)
    is opened once and never re-read.
    """

    def __init__(self):
        self._stack = None
        self._zf = None
        self._key = None
        self._dirs = set()

    def _open(self, rec):
        key = (str(rec["zip_path"]), "/".join(rec["nested"]))
        if key == self._key and self._zf is not None:
            return
        self.close()
        stack = contextlib.ExitStack()
        try:
            zf = stack.enter_context(zipfile.ZipFile(rec["zip_path"], "r"))
            # Descend the nested chain, spooling each inner archive to disk.
            for component in rec["nested"]:
                fh = stack.enter_context(_spool_member(zf, zf.getinfo(component)))
                zf = stack.enter_context(zipfile.ZipFile(fh))
        except Exception:
            stack.close()
            raise
        self._stack, self._zf, self._key = stack, zf, key

    def ensure_dir(self, directory: Path):
        key = str(directory)
        if key not in self._dirs:
            Path(directory).mkdir(parents=True, exist_ok=True)
            self._dirs.add(key)

    def copy_to(self, rec, dest: Path):
        """Byte-for-byte copy of a member (written via .part, then renamed)."""
        self._open(rec)
        self.ensure_dir(dest.parent)
        tmp = dest.with_name(dest.name + ".part")
        try:
            with self._zf.open(rec["member"]) as src, open(tmp, "wb") as out:
                shutil.copyfileobj(src, out, COPY_CHUNK)
                out.flush()
            os.replace(tmp, dest)
        finally:
            tmp.unlink(missing_ok=True)

    def spool(self, rec):
        """Seekable copy of a member (RAM up to SPOOL_LIMIT, then a temp file).

        Needed when the member has to be parsed and resaved; caller closes it.
        """
        self._open(rec)
        tmp = tempfile.SpooledTemporaryFile(max_size=SPOOL_LIMIT)
        try:
            with self._zf.open(rec["member"]) as src:
                shutil.copyfileobj(src, tmp, COPY_CHUNK)
            tmp.seek(0)
        except Exception:
            tmp.close()
            raise
        return tmp

    def close(self):
        if self._stack is not None:
            self._stack.close()
        self._stack = self._zf = self._key = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _patient_id_of(fp):
    """(0010,0020) held by a seekable DICOM stream: value, '' if absent, None if unreadable."""
    fp.seek(0)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ds = pydicom.dcmread(fp, force=True, stop_before_pixels=True,
                                 specific_tags=["PatientID"])
    except Exception:
        return None
    return str(getattr(ds, "PatientID", "")).strip()


def _write_fixed_patient_id(fp, dest: Path, patient_id: str):
    """Resave the member in `fp` to `dest` with (0010,0020) PatientID = `patient_id`.

    Nothing else is changed: no other tag is added, edited or removed, the SOP
    Instance UID, the SOP Class UID and the transfer syntax are left exactly as
    read. The dataset is read in full (pixel data included) because `save_as`
    needs a complete dataset to write - one file in RAM at a time, and the file
    is only put in place after it has been read back and checked.
    """
    fp.seek(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = pydicom.dcmread(fp, force=True)
    ds.PatientID = patient_id

    tmp = dest.with_name(dest.name + ".part")
    problems = []
    for kwargs in (
        {"enforce_file_format": True},                      # conformant Part-10 output
        {"enforce_file_format": False},                     # source had no meta/preamble
        {"enforce_file_format": True, "force_encoding": True},  # odd SpecificCharacterSet
    ):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with open(tmp, "wb") as out:
                    ds.save_as(out, **kwargs)
            with open(tmp, "rb") as fh:                     # verify before publishing
                got = _patient_id_of(fh)
            if got != patient_id:
                raise ValueError(f"PatientID read back as {got!r}")
            os.replace(tmp, dest)
            return
        except Exception as exc:
            problems.append(f"{exc}")
    tmp.unlink(missing_ok=True)
    raise RuntimeError("could not write corrected PatientID (" + "; ".join(problems[:2]) + ")")


def _origin_key(rec):
    """Stable identity of the physical source of a record."""
    return (str(rec["zip_path"]), "/".join(rec["nested"]), rec["member"])


def _instance_key(rec):
    """Identity of the instance on disk: the SOP Instance UID is unique.

    Two archives/members carrying the same patient + study + series + SOP UID are
    the same instance exported twice (verified on this corpus: the only differing
    element is (0008,0013) InstanceCreationTime), so one file must win.
    """
    return (rec["cleaned_id"], rec["study_uid"], rec["series_uid"], rec["sop_uid"])


def unique_instances(records):
    """Collapse duplicate copies of the same SOP instance, last copy in archive
    order wins (i.e. it overwrites the earlier one) - deterministic per run."""
    winners = {}
    for rec in sorted(records, key=_origin_key):
        winners[_instance_key(rec)] = rec
    return sorted(winners.values(), key=_origin_key)


def _joined(values):
    return "; ".join(sorted({str(v) for v in values if v}))


def _archive_label(rec):
    if rec["nested"]:
        return f"{rec['zip_path'].name}!{'>'.join(rec['nested'])}"
    return rec["zip_path"].name


def _instance_path(rec, out_root: Path) -> Path:
    return Path(out_root) / rec["cleaned_id"] / rec["study_uid"] / rec["series_uid"] / rec["instance_name"]


def write_records(records, out_root: Path, reader: MemberReader, label="files", progress=True):
    """Write records to <out_root>/<pid>/<study>/<series>/<sop>_<modality>.dcm.

    One file per SOP instance: duplicate copies of the same instance are skipped
    (the surviving copy is the one that would have overwritten the others).

    (0010,0020) PatientID is corrected to the cleaned id the files were matched
    on; files already carrying that id are copied byte for byte, untouched.
    """
    empty = {"written": 0, "bytes": 0, "duplicates": 0, "failed": 0,
             "patched": 0, "per_archive": {}}
    if not records:
        return empty

    # Archive order: each zip is opened exactly once, sequential I/O only.
    originals = len(records)
    records = unique_instances(records)
    duplicates = originals - len(records)
    every = max(500, len(records) // 40)   # ~40 progress lines at most
    started = time.time()
    written = total_bytes = failed = 0
    patched = verbatim = 0
    per_archive = {}

    for i, rec in enumerate(records, 1):
        dest = _instance_path(rec, out_root)
        try:
            if rec["patient_id"] == rec["cleaned_id"] and rec["patient_id"]:
                reader.copy_to(rec, dest)          # id already correct: raw bytes
                verbatim += 1
            else:
                # the file's PatientID is not the id it was matched with (it had
                # alphabets / separators, or was missing): write the corrected one
                reader.ensure_dir(dest.parent)
                with reader.spool(rec) as fp:
                    _write_fixed_patient_id(fp, dest, rec["cleaned_id"])
                patched += 1
        except Exception as exc:  # a single bad member must not kill the run
            print(f"  ! could not write {rec['source']}: {exc}")
            failed += 1
            continue
        written += 1
        total_bytes += rec["size"]
        stats = per_archive.setdefault(_archive_label(rec), [0, 0])
        stats[0] += 1
        stats[1] += rec["size"]
        if progress and i % every == 0:
            print(f"    {label}: {i}/{len(records)} ({100 * i // len(records)}%, {time.time() - started:.0f}s)", flush=True)
    if progress and len(records) % every:
        print(f"    {label}: {len(records)}/{len(records)} (100%, {time.time() - started:.0f}s)", flush=True)

    print(f"    {label}: {written} files, {total_bytes / 1e9:.2f} GB"
          f"{f', {duplicates} duplicate SOP instances skipped' if duplicates else ''}"
          f"{f', {failed} failed' if failed else ''} in {time.time() - started:.0f}s")
    if written:
        print(f"      PatientID (0010,0020): {patched} corrected to the matched id,"
              f" {verbatim} already correct")
    for name, (count, nbytes) in sorted(per_archive.items()):
        print(f"      {name}: {count} files, {nbytes / 1e9:.2f} GB")
    return {"written": written, "bytes": total_bytes, "duplicates": duplicates,
            "failed": failed, "patched": patched, "per_archive": per_archive}


# --------------------------------------------------------------------------
# CSV handling
# --------------------------------------------------------------------------
def process_csv_file(src: Path, id_column: str, cleaned_dst: Path):
    """Copy the CSV to `cleaned_dst` with a cleaned patient id column.

    Returns a list of cleaned patient records:
        {"original_id": str, "cleaned_id": str}
    """
    cleaned_records = []

    with open(src, "r", newline="", encoding="utf-8-sig") as fin, \
         open(cleaned_dst, "w", newline="", encoding="utf-8") as fout:

        reader = csv.DictReader(fin)
        fieldnames = list(reader.fieldnames or [])
        writer = csv.writer(fout)
        writer.writerow(fieldnames)

        for row in reader:
            values = {}
            for key in fieldnames:
                val = row.get(key)
                values[key] = "" if val is None else val

            original_id = str(values.get(id_column, "") or "").strip()
            cleaned_id = normalize_id(original_id)

            values[id_column] = cleaned_id
            writer.writerow([values[k] for k in fieldnames])

            cleaned_records.append({
                "original_id": original_id,
                "cleaned_id": cleaned_id,
            })

    return cleaned_records


def select_id_column(fieldnames):
    """Prompt the user to pick the CSV column that holds the patient id."""
    print("\nAvailable columns in the CSV:")
    for i, col in enumerate(fieldnames):
        print(f"  [{i}] {col}")
    while True:
        try:
            choice = input("Select the column that corresponds to the patient id (enter number or column name): ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n  No selection made - skipping this CSV.")
            return None
        if choice in fieldnames:
            return choice
        try:
            idx = int(choice)
            if 0 <= idx < len(fieldnames):
                return fieldnames[idx]
        except ValueError:
            pass
        print(f"Invalid selection: '{choice}'. Please try again.")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def collect_csv_jobs(csv_files, cleaned_dir):
    """Prompt for the id column and clean every CSV.

    Done *before* the DICOM scan so the script is interactive (and visibly
    alive) the moment it starts, instead of looking hung while it loads the
    whole DICOM corpus.
    """
    jobs = []
    for csv_path in csv_files:
        print(f"\nProcessing CSV file: {csv_path}")
        with open(csv_path, "r", newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            print("  CSV appears to have no header row - skipping.")
            continue

        id_column = select_id_column(fieldnames)
        if id_column is None:
            continue
        print(f"  Using column '{id_column}' as patient id column.")

        if id_column not in fieldnames:
            print(f"  Error: selected column '{id_column}' was not found in the CSV header. Skipping.")
            continue

        cleaned_dst = cleaned_dir / csv_path.name
        cleaned_records = process_csv_file(csv_path, id_column, cleaned_dst)
        print(f"  Cleaned CSV copy written to: {cleaned_dst}")

        # Unique cleaned patient ids from this CSV (order preserved, no duplicates)
        unique_pids = []
        originals = {}
        for rec in cleaned_records:
            pid = rec["cleaned_id"]
            if pid not in originals:
                originals[pid] = rec["original_id"]
                unique_pids.append(pid)

        jobs.append({
            "csv_path": csv_path,
            "unique_pids": unique_pids,
            "originals": originals,
        })
    return jobs


def report_csv_job(job, dicom_index, results_dir, matched_dir, output_dir):
    """Write <stem>_matched.csv / <stem>_not_found.csv for one input CSV.

    Returns the flat list of DICOM records that still have to be copied.
    Rows are streamed straight to disk - nothing accumulates in memory.
    """
    csv_path = job["csv_path"]
    originals = job["originals"]
    stem = csv_path.stem

    matched_path = results_dir / f"{stem}_matched.csv"
    to_copy = []
    n_matched = 0
    with open(matched_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["PatientId", "StudyInstanceUIDs", "SeriesInstanceUIDs",
                        "Modalities", "Instances", "FolderPaths"])
        for pid in job["unique_pids"]:
            recs = dicom_index.get(pid)
            if recs is None:
                continue
            # ---- MATCHED (metadata only, copying happens in pass 2) ----
            series = {}
            for r in recs:
                series.setdefault((r["study_uid"], r["series_uid"]), True)

            folder_paths = sorted(
                str((Path(matched_dir) / pid / study / sid).relative_to(output_dir))
                for study, sid in series
            )
            writer.writerow([
                pid,
                _joined(s for s, _ in series),
                _joined(s for _, s in series),
                _joined(r["modality"] for r in recs),
                len(unique_instances(recs)),
                "; ".join(folder_paths),
            ])
            to_copy.extend(recs)
            n_matched += 1

    not_found_path = results_dir / f"{stem}_not_found.csv"
    n_not_found = 0
    with open(not_found_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["PatientId", "OriginalId"])
        for pid in job["unique_pids"]:
            if pid in dicom_index:
                continue
            writer.writerow([pid, originals.get(pid, "")])
            n_not_found += 1

    print(f"  Matched patients: {n_matched}")
    print(f"  Not found patients: {n_not_found}")
    print(f"  Matched results CSV: {matched_path}")
    print(f"  Not found results CSV: {not_found_path}")
    return to_copy


def write_unmatched_report(dicom_index, all_seen, results_dir):
    """Write results/unmatched.csv: which DICOM patients matched no CSV row.

    Counts are unique SOP instances (what actually lands in unmatched/), and the
    source column names the archive (and inner archive, if nested) they came from.
    """
    path = Path(results_dir) / "unmatched.csv"
    patients = [(pid, recs) for pid, recs in sorted(dicom_index.items()) if pid not in all_seen]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["PatientId", "Instances", "Modalities", "SourceArchives"])
        for pid, recs in patients:
            writer.writerow([
                pid,
                len(unique_instances(recs)),
                _joined(r["modality"] for r in recs),
                _joined(_archive_label(r) for r in recs),
            ])
    return path, len(patients)


def main():
    parser = argparse.ArgumentParser(description="Match DICOM patient ids to cleaned CSV patient ids.")
    parser.add_argument("--csv-dir", default=DEFAULT_CSV_DIR, help=f"Folder with input CSV files (default: {DEFAULT_CSV_DIR})")
    parser.add_argument("--dicom-dir", default=DEFAULT_DICOM_DIR, help=f"Folder with DICOM zip files (default: {DEFAULT_DICOM_DIR})")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help=f"Folder for cleaned / matched / unmatched outputs (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--skip-unmatched", action="store_true",
                        help="Do not copy unmatched DICOM files (saves a full second pass and one copy of the corpus)")
    parser.add_argument("--no-progress", action="store_true", help="Quiet mode: only print per-file summaries")
    parser.add_argument("--force", action="store_true", help="Continue even if the output drive looks too small")
    args = parser.parse_args()

    progress = not args.no_progress

    csv_dir = Path(args.csv_dir)
    dicom_dir = Path(args.dicom_dir)
    output_dir = Path(args.output_dir)

    cleaned_dir = output_dir / "cleaned"
    matched_dir = output_dir / "matched"
    unmatched_dir = output_dir / "unmatched"
    results_dir = output_dir / "results"

    # Start from a clean output state on every run
    for sub in (cleaned_dir, matched_dir, unmatched_dir, results_dir):
        if sub.exists():
            print(f"Clearing {sub} ...")
            shutil.rmtree(sub)
        sub.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(csv_dir.glob("*.csv"))
    if not csv_files:
        print(f"No CSV files found in '{csv_dir}'.")
        return

    # ---- 1. CSVs first: prompt the user, write cleaned copies ----
    jobs = collect_csv_jobs(csv_files, cleaned_dir)
    if not jobs:
        print("No usable CSV files - nothing to do.")
        return

    # ---- 2. Pass 1: DICOM tags only (metadata index, tiny memory footprint) ----
    dicom_index = build_dicom_index(dicom_dir, progress=progress)

    # ---- 3. Reports per CSV, from metadata only ----
    to_copy = []
    seen_ids = set()
    for job in jobs:
        for rec in report_csv_job(job, dicom_index, results_dir, matched_dir, output_dir):
            key = _origin_key(rec)
            if key not in seen_ids:      # same patient listed in several CSVs
                seen_ids.add(key)
                to_copy.append(rec)

    all_seen = {pid for job in jobs for pid in job["unique_pids"]}
    # Unmatched = DICOM patients that appear in no CSV at all (one pass, not one
    # per CSV file as before).
    unmatched = [r for pid, recs in sorted(dicom_index.items()) if pid not in all_seen for r in recs]
    unmatched_report, n_unmatched_patients = write_unmatched_report(dicom_index, all_seen, results_dir)

    # One file per SOP instance: duplicated SOP UIDs are the same instance exported
    # twice, so the duplicate copies are dropped instead of written twice.
    matched_records = unique_instances(to_copy)
    unmatched_records = unique_instances(unmatched)
    dup_matched = len(to_copy) - len(matched_records)
    dup_unmatched = len(unmatched) - len(unmatched_records)
    print(f"\n  Unmatched patients: {n_unmatched_patients} ({len(unmatched_records)} DICOM files)")
    print(f"  Unmatched results CSV: {unmatched_report}")
    if dup_matched or dup_unmatched:
        print(f"  Duplicate SOP instances (same instance exported twice, one file kept):"
              f" {dup_matched} matched + {dup_unmatched} unmatched copies skipped")

    # ---- 4. Disk budget check before touching any bytes ----
    needed = sum(r["size"] for r in matched_records)
    if not args.skip_unmatched:
        needed += sum(r["size"] for r in unmatched_records)
    free = shutil.disk_usage(output_dir).free
    print(f"\nNeed {needed / 1e9:.1f} GB of output space, {free / 1e9:.1f} GB free on '{output_dir}'")
    if needed > free and not args.force:
        print("Aborting before writing: not enough free space (use --force to override).")
        return

    # ---- 5. Pass 2: stream the files out, one archive open at a time ----
    print(f"Pass 2/2: writing {len(matched_records)} matched and {len(unmatched_records)} unmatched DICOM files")
    with MemberReader() as reader:
        write_records(matched_records, matched_dir, reader, "matched", progress=progress)
        if not args.skip_unmatched:
            write_records(unmatched_records, unmatched_dir, reader, "unmatched", progress=progress)
        else:
            print("  Skipping unmatched DICOM copy (--skip-unmatched)")

    print("\nDone.")
    print(f"  Matched DICOM files: {matched_dir}")
    if not args.skip_unmatched:
        print(f"  Unmatched DICOM files: {unmatched_dir}")
    print(f"  Result CSVs: {results_dir} (matched / not_found per CSV + unmatched.csv)")


if __name__ == "__main__":
    main()
