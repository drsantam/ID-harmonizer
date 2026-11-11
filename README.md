# ID Haromonizer

A CLI tool that matches patient ids from CSV files to DICOM archives and
copies the matching DICOM files into a clean, organized output folder.

## What it does

1. Reads CSV file(s) from the `csv/` folder.
2. Prompts you to select the column that holds the patient id.
3. Writes a cleaned copy of the CSV to the output folder (patient id stripped
   of alphabets, special characters/spaces converted to underscores).
4. Reads all DICOM files inside zip files (nested zips supported) AND all loose
   DICOM files in subfolders of the `dicom/` folder (scanned recursively),
   extracting the patient ids.
5. Matches DICOM ids to the cleaned CSV ids and copies matching files to
   `output/matched/`, organized as `Study > Series > Instance`, named
   `<SOP Instance UID>_<Modality>.dcm`. The `PatientID` tag of each copied
   file is set to the cleaned id it was matched on.
6. Keeps unmatched DICOM files in `output/unmatched/`.
7. Writes result CSVs (not-found patients, matched patients, unmatched
   DICOM patients) to `output/results/`.

Original input CSV files are never modified.

## Requirements

- Python 3
- Dependencies (see `requirements.txt`): `pydicom`, `pandas`, `numpy`

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Running

Place your input CSV files in `csv/` and your DICOM inputs in `dicom/` (zip
files, plus subfolders containing loose DICOM files - subfolders are scanned
recursively), then run:

```bash
python3 id_haromonizer.py
```

You will be prompted to select the patient-id column (enter the column number
or name). Progress is printed continuously throughout the run.

### Options

```
--csv-dir DIR      folder with input CSV files (default: csv)
--dicom-dir DIR    folder with DICOM zip files (default: dicom)
--output-dir DIR   folder for cleaned / matched / unmatched outputs (default: output)
--skip-unmatched   don't copy unmatched DICOM files to the output folder
--no-progress      quiet mode: only print per-file summaries
--force            continue even if the output drive looks too small
```

Example with custom directories:

```bash
python3 id_haromonizer.py --csv-dir /data/csv --dicom-dir /data/dicom --output-dir /data/out
```
