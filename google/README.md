# google

## contains

_in ./gcp.py_

- mvFiles: move files to a folder.
- lsFiles: list all files.
- cpFiles: copy many files to a foler.
- catFiles: get data in many files.
- rmFiles: remove many files.
- recoverFiles: if bucket has versioning enabled, retrieve list of files that have been deleted.
- patternRN: following a renaminig dict, rename a bunch of files in a set of locations.
- get_all_sizes: get file sizes in a folder.
- exists: given a list of file paths, get if files exist or not.
- extractSize: extract file size from ls command.
- extractPath: extract file path from ls command.
- extractHash: extract file hash from ls command.

_in google\_sheet.py_

- dfToSheet: uploads a given dataframe to a given googleSheet location

GSheet (class) *WIP*
  - get_last_modified_date
  - get_size
  - read_sheet
  - read_row
  - read_column
  - write_column

  # highly recommended:

  - gsutil
  - pygsheet