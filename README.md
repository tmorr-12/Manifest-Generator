# Manifest Generator

Wellcome to Manifest Generator, the swiss army knife of manifest generation! Manifest generator is a command-line tool for generating sample manifests from FASTQ files, supporting short, long, and hybrid sequencing read types.

---

## Requirements

- Python 3.10+
- `pandas 2.3+`
- [`lrge`](https://github.com/mbhall88/lrge) _(optional, required for genome size estimation)_

---

## Install

> **Warning:** `lrge` is not installed via pip. To install it, please see the [`lrge` repo](https://github.com/mbhall88/lrge) linked above.

```bash
python3 -m venv .venv \
source .venv/bin/activate \
pip install requirements.txt
```

---

## Usage

```bash
manifest_generator.py [--from_dir | --from_dir_recursive | --from_paths | --from_paths_id] [flags]
```

---

## Modes

### `short`

The short mode is used by default. Pairs R1 and R2 short read files by sample ID. Sample ID is inferred from filename (see [here](#id-parsing))

Output columns: `ID`, `R1`, `R2`

### `long`

Assigns long reads to sample IDs.

Output columns: `ID`, `long_fastq`, `genome_size`

### `hybrid`

Combines short read pairs with a long read per sample.

Output columns: `ID`, `R1`, `R2`, `long_fastq`, `genome_size`

---

## Input Sources

| Source                                   | Description                                                                             |
| ---------------------------------------- | --------------------------------------------------------------------------------------- |
| `--from_dir <dir> [<dir> ...]`           | Collect FASTQs by searching one or more directories                                     |
| `--from_dir_recursive <dir> [<dir> ...]` | Collect FASTQs by recursively searching one or more directories(requires `--max_depth`) |
| `--from_paths <file>`                    | Supply a list of FASTQ file paths from which to build a manifest                        |
| `--from_paths_id <file>`                 | Supply a csv containing paired read paths with ID                                       |

---

## Flags

| Flag                   | Default        | Description                                             |
| ---------------------- | -------------- | ------------------------------------------------------- |
| `-m`, `--mode`         | `short`        | Manifest type: `short`, `long`, or `hybrid`             |
| `-o`, `--outdir`       | `./`           | Output directory                                        |
| `-n`, `--name`         | `manifest.csv` | Now YOU can name your very own manifest!                |
| `-d`, `--max_depth`    | `None`         | Max directory depth for `--from_dir_recursive`          |
| `--genome_size`        | `False`        | Estimate genome size for long reads using `lrge`        |
| `-t`, `--threads`      | `1`            | Threads for genome size estimation                      |
| `--duplicate_handling` | `error`        | How to handle duplicate IDs: `error`, `warn`, or `drop` |
| `-v`, `--verbose`      | `False`        | Enable debug logging                                    |

---

## Duplicate Handling

The `--duplicate_handling` flag controls behaviour when duplicates are found:

| Method  | Behaviour                                     |
| ------- | --------------------------------------------- |
| `error` | Log error and exit (default)                  |
| `warn`  | Log warning and continue                      |
| `drop`  | Drop duplicates, keep first occurrence per ID |

Recommended usage:

| Method  | Behaviour                                                                             |
| ------- | ------------------------------------------------------------------------------------- |
| `error` | When you are unsure if duplicates are present                                         |
| `warn`  | When you are unsure if duplicates are present and wish to inspect the output yourself |
| `drop`  | When you don't care which read is being used                                          |

---

## ID parsing

IDs will only be parsed correctly from filenames with the following rules:

- Filenames must match one of the following suffixes: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`
- For the `short` mode, presence of `_R1` or `_1` are used to infer R1 files, while `_R2` or `_2` are used to infer R2 files.
  > :warning: IDs will not be correctly extracted from filenames that contain these patterns (`_R?[12]`) internally, e.g. `file_1_1.fastq.gz`.

---

## Examples

```bash
# Short read manifest from a directory
manifest_generator.py --from_dir /path/to/fastqs -m short -o ./output

# Hybrid manifest, searching recursively up to depth 3
manifest_generator.py --from_dir_recursive /path/to/project -d 3 -m short -o ./output

# Long read manifest with genome size estimation
manifest_generator.py --from_dir /path/to/fastqs -m long --genome_size -t 8

# From a list of paths, warn on duplicates
manifest_generator.py --from_paths reads.txt --duplicate_handling warn

# Hybrid read manifest
manifest_generator.py --from_paths_id paths_id.csv -m hybrid -o ./output
```

---

## Output

A `.csv` manifest is written to `--outdir` with the filename specified by `--name`. A timestamped log file is also written to the same directory.

```
output/
├── manifest.csv
└── 2025-01-01_12-00-00.log
```
