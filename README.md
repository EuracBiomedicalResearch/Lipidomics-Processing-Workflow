# Lipidomics-Processing_Workflow

**Semi-automated workflow for High-confidence lipid annotation using a public SRM 1950-derived lipid database**

A vendor-independent, open-source workflow for Untargeted LC-MS-based lipidomics, data processing and annotation.

---

## Overview

This workflow provides:

- **Automated RT adjustment** using a Lipid reference set (containing SPLASH LIPIDOMIX and some endogenous compounds)
- **Lipid annotation** against a curated SRM1950 human plasma lipid database
- **Isotope pattern validation** for high-confidence identification
- **Adduct profile matching** for improved annotation
- **Internal standard normalization** by lipid subclass
- **Quality control filtering** based on QC sample RSD

### Supported Data

- **Ionization modes**: Positive (POS) and Negative (NEG)
- **Instrument platforms**: Vendor-independent (requires .mzML format)
- **Sample types**: QC, M-VAMS, W-DBS, C-qDBS, Plasma (configurable)

---

## Project Structure

Each analysis is split by **ionization mode** — positive and negative each
have their own Quarto documents so that parameter choices are fully traceable.

```
CEMBIO-EURAC/
│
├── generic_workflow/                # Reusable workflow templates
│   ├── LipidDatabase_R.xlsx         #   SRM 1950 Lipid database (sheet 4=POS, sheet 5=NEG)
│   ├── positive/                    #   Positive ionization mode templates
│   │   ├── Preprocessing_pos.qmd    #   Stage 1 (POLARITY = "pos")
│   │   └── Annotation_pos.qmd       #   Stage 2 (POLARITY = "pos")
│   ├── negative/                    #   Negative ionisation mode templates
│   │   ├── Preprocessing_neg.qmd    #   Stage 1 (POLARITY = "neg")
│   │   └── Annotation_neg.qmd       #   Stage 2 (POLARITY = "neg")
│   └── POS_NEG_merge.qmd            #   Stage 3: Cross-ionization mode Data Integration
│
├── applications/                    # Study-specific applications (self-contained)
│   ├── MICROSAMPLING_study/
│   │   ├── LipidDatabase_R.xlsx     #   SRM 1950 Lipid database (shared by POS + NEG)
│   │   ├── positive/                #   Positive ionization mode analysis
│   │   │   ├── Preprocessing_pos.qmd   #     Configured for MICROSAMPLING, POS
│   │   │   ├── Annotation_pos.qmd      #     Configured for MICROSAMPLING, POS
│   │   │   ├── seq_pos_MICROSAMPLING.xlsx   #     Sample sequence
│   │   │   ├── pos_lipid_reference_set.xlsx  # Reference lipids
│   │   │   └── data/                #     .mzML files (gitignored)
│   │   ├── negative/                #   Negative ionization mode analysis
│   │   │   ├── Preprocessing_neg.qmd   #     Configured for MICROSAMPLING, NEG
│   │   │   ├── Annotation_neg.qmd      #     Configured for MICROSAMPLING, NEG
│   │   │   ├── seq_neg_MICROSAMPLING.xlsx   #     Sample sequence
│   │   │   ├── neg_lipid_reference_set.xlsx  # Reference lipids
│   │   │   └── data/               #     .mzML files (gitignored)
│   │   └── POS_NEG_merge.qmd       #   Stage 3: Cross-ionization mode Data Integration MICROSAMPLING study
│   └── METFORMIN-HIIE_study/        #   (same structure as MICROSAMPLING study)
│       ├── positive/ ...
│       ├── negative/ ...
│       └── POS_NEG_merge.qmd
│
├── R/                               # Shared helper functions
│   ├── lipid_helpers.R              #   All reusable functions
│   └── create_sqlite_database.R     #   SQLite database creation utility
│
└── README.md                        # This file
```

---

## Quick Start

### Docker setup and usage

Docker bundles R, Quarto, compilers and the locked R packages. With Docker
installed and running, you can use this option instead of the native Conda/R
setup below. Run the commands from the repository root.

Build the image and check its environment:

```bash
docker build --platform linux/amd64 -t cembio-eurac:main-docker .
docker run --rm --platform linux/amd64 cembio-eurac:main-docker check
```

The first build downloads and compiles the dependencies and can take a long
time. It also checks native compilation and renders a Quarto smoke test.
The image targets Linux x86-64; Apple Silicon uses emulation through the
`--platform linux/amd64` option and can be slower.

MICROSAMPLING is configured to download public MetaboLights data. For METFORMIN,
provide `METFORMIN-HIIE_pos.sqlite` in
`applications/METFORMIN-HIIE_study/positive/data/` and
`METFORMIN-HIIE_neg.sqlite` in its `negative/data/` folder. If these databases
are absent, the workflow tries to create them from the `.mzML` files named in
the corresponding sequence workbook; place those files in each mode's `data/`
folder. The other required inputs are described in **Prepare Input Files** below.

Run a complete study with the checkout mounted at `/workspace` and a persistent
download cache:

```bash
docker run --rm --init --platform linux/amd64 \
  --mount "type=bind,source=$(pwd),target=/workspace" \
  --mount type=volume,source=cembio-cache,target=/cache \
  cembio-eurac:main-docker MICROSAMPLING

docker run --rm --init --platform linux/amd64 \
  --mount "type=bind,source=$(pwd),target=/workspace" \
  --mount type=volume,source=cembio-cache,target=/cache \
  cembio-eurac:main-docker METFORMIN
```

On Linux, add `--user "$(id -u):$(id -g)"` before the image name, starting with
your first run, so generated files belong to your account. Keep the same user,
cache volume and `/workspace` mount path across stages. Docker Desktop users
should ensure the checkout is shared with Docker. These examples use Bash
syntax, including WSL on Windows.

Results, HTML reports, SQLite files and annotation workbooks are written into
the mounted checkout and remain after the container exits. Downloaded data
remain in the `cembio-cache` volume. Raw data, previous results and local R
libraries are excluded from the image. Saved experiments can refer to backing
files by their container paths, so create and reuse preprocessing outputs with
the same mounts.

To run individual stages, replace the study arguments at the end of the
`docker run` command:

| Arguments | Action |
|-----------|--------|
| `MICROSAMPLING preprocessing` | Preprocess both polarities |
| `MICROSAMPLING preprocessing pos` | Positive preprocessing only |
| `MICROSAMPLING preprocessing neg` | Negative preprocessing only |
| `MICROSAMPLING annotation pos` | Annotate saved positive preprocessing |
| `MICROSAMPLING annotation neg` | Annotate saved negative preprocessing |
| `MICROSAMPLING merge` | Merge saved annotations from both polarities |
| `help` | Show available commands |

Both study names support the same stages. Without a stage argument, the runner
executes positive preprocessing and annotation, then negative preprocessing
and annotation, then the merge. Annotation requires preprocessing outputs;
merge requires annotation outputs from both modes. Rerunning a stage overwrites
its generated outputs.

Edit scientific parameters, including `CORES_NB`, in the study's `.qmd` files
before running. Use paths under `/workspace` for custom inputs. Allocate Docker
enough CPU, memory and disk for the LC-MS dataset. Mounted workflow edits take
effect on the next run; rebuild the image after changing `environment.yml` or
`renv.lock`.

The [validation checkpoints](#user-validation-checkpoints) need scientific
review; an unattended render does not pause at them. For manual curation, save
reviewed ambiguity workbooks under separate filenames and update the annotation
document's curation reads to use those copies: its export chunk regenerates the
default ambiguity workbooks on each render.

### 1. The Reproducible Environment

The project uses two complementary environment layers:

- Conda supplies R 4.6.0, Quarto, compilers, and native libraries from
  `environment.yml`.
- `renv` supplies the exact CRAN, Bioconductor 3.23, and GitHub package
  versions recorded in `renv.lock`.

From the repository root, create the Conda environment:

```bash
mamba env create --file environment.yml
```

If Mamba is unavailable, use `conda env create --file environment.yml`
instead. This step is only required **once**.

Activate the Conda environment. This step must be done each time **before**
running R or Quarto. Eventually deactivate the conda environment after running
the workflows to restore the default system setup (using `conda deactivate`).

```bash
conda activate cembio_eurac
```

Set up or restore the R package library and verify the complete installation:

```bash
Rscript scripts/bootstrap_environment.R
Rscript scripts/check_environment.R
```

The first time these commands are executed, all required R packages are
installed. Any subsequent call will restore the cached libraries.

Starting R in the base folder of the repository will automatically set up and
load the environment (pre-configured by the *renv.lock* file).

### 2. Prepare Input Files

Each ionization mode folder is self-contained. Place these files inside your
polarity folder (e.g. `applications/my_study/positive/`):

1. **LC-MS raw data**: Place `.mzML` files in the `data/` subfolder

2. **Sample sequence file** (`seq_pos_<study_id>.xlsx` or `seq_neg_<study_id>.xlsx`):

   | file_name | sample_name | sample_type | injection_index |
   |-----------|-------------|-------------|------------------|
   | D01P_pos.mzML | D01P | Plasma | 1 |
   | QC_1_pos.mzML | QC_1 | QC | 2 |

3. **SRM 1950-derived Lipid database**: `LipidDatabase_R.xlsx` in the study folder (one level above each polarity folder)

4. **Lipid Reference Set**: `pos_lipid_reference_set.xlsx` or `neg_lipid_reference_set.xlsx`

### 3. Run the Workflow

Each polarity has dedicated Quarto documents with hardcoded `POLARITY`,
so parameter choices are fully traceable.

| Step | Files | Description |
|------|-------|-------------|
| 1 | `positive/Preprocessing_pos.qmd`, `negative/Preprocessing_neg.qmd` | Data import and validation, peak detection, RT alignment, correspondence and gap filling |
| 2 | `positive/Annotation_pos.qmd`, `negative/Annotation_neg.qmd` | SRM1950 Database RT adjustment, Multi-evidence annotation, normalization, QC |
| 3 | `POS_NEG_merge.qmd` | Positive and negative ionization mode integration for coverage reporting and downstream analysis |

The MICROSAMPLING preprocessing documents retrieve their mzML files from the
public MetaboLights study `MTBLS10722`. The first positive and negative
preprocessing runs download 48 files per polarity (approximately 2 GB total)
into the user's BiocFileCache, normally `~/.cache/R/BiocFileCache`. Later runs
reuse the cached files.

To ensure the reproducible R environment setup from step 1 is used, start R in
the base directory of the repository.

**Example: run the MICROSAMPLING study**

1. Activate the Conda environment set up in step 1: `conda activate cembio_eurac`.
2. Start R in the base directory of the repository and render the documents in
   this order:

   ```r
   quarto::quarto_render("applications/MICROSAMPLING_study/positive/Preprocessing_pos.qmd")
   quarto::quarto_render("applications/MICROSAMPLING_study/positive/Annotation_pos.qmd")
   quarto::quarto_render("applications/MICROSAMPLING_study/negative/Preprocessing_neg.qmd")
   quarto::quarto_render("applications/MICROSAMPLING_study/negative/Annotation_neg.qmd")
   quarto::quarto_render("applications/MICROSAMPLING_study/POS_NEG_merge.qmd")
   ```

### 4. Start a New Study

1. Create a new folder under `applications/` (e.g. `applications/my_study/`)
2. Create `positive/` and `negative/` subfolders
3. Copy templates from `generic_workflow/positive/` and `generic_workflow/negative/`
4. Copy `generic_workflow/POS_NEG_merge.qmd` to the study root
5. In each file, set:
   - `PROJECT_ROOT <- "../../.."` (path back to project root)
   - `STUDY_ID <- "my_study"`
   - Adjust study-specific parameters (sample types, colors, volume factors, etc.)
6. Copy `LipidDatabase_R.xlsx` into the study folder (one level above `positive/`/`negative/`)
7. Add polarity-specific files to each polarity folder:
   - `seq_pos_my_study.xlsx` / `seq_neg_my_study.xlsx` (sample sequences)
   - `pos_lipid_reference_set.xlsx` / `neg_lipid_reference_set.xlsx` (copy from another study)
   - `data/` subfolder with `.mzML` files

---

## User Validation Checkpoints

The workflow includes several interactive checkpoints:

| Checkpoint | Location | Action Required |
|------------|----------|-----------------|
| RT Filter Range | Preprocessing | Adjust RT filter based on BPC |
| Lipid Reference Set EIC | Preprocessing | Verify internal standard signals |
| Peak Detection Results | Preprocessing | Check peak detection quality |
| RT Adjustment Model | Annotation | Validate R-squared and residuals |
| Isotope Pattern Validation | Annotation | Review mirror plots |
| Manual Ambiguity Resolution | Annotation | Curate Excel files |

---

## Configuration Parameters

Key parameters in the Configuration section of each workflow file:

```r
PROJECT_ROOT <- "../../.."  # Path to project root (for shared R/ helper functions)
STUDY_ID <- "MICROSAMPLING" # Study prefix for saved objects and seq file names
POLARITY <- "pos"           # Hardcoded per folder ("pos" in positive/, "neg" in negative/)
CORES_NB <- 4               # Parallel processing cores
PPM <- 10                   # m/z tolerance for peak detection
MATCH_PPM <- 20             # m/z tolerance for database matching
MATCH_RT_TOL <- 20          # RT tolerance for matching (seconds)
ISOPEAK_SIM_THRESHOLD <- 0.78  # Minimum isotope similarity
RSD_THRESHOLD <- 0.3        # QC RSD filter (30%)
```

---

## 🆘 Troubleshooting

- The `quarto_render()` call does not run R from the configured environment
  (*Quick start*, point 1.): check if an environment variable `QUARTO_R` is set
  (e.g. using `Sys.getenv("QUARTO_R")` in R or `echo $QUARTO_R` in a shell) and
  if so, *unset* it.

## References

1. **SRM 1950-DERIVED Lipid Database**: <https://doi.org/10.1016/j.jlr.2024.100671>
2. **Original study (MICROSAMPLING dataset)**: <https://doi.org/10.1016/j.microc.2025.113760>
3. **Application study (METFORMIN-HIIE dataset)**: <https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00480>

---

## Contributors

- CEMBIO-EURAC Team
- Sara Londoño-Osorio
- Philippe Louail

---

## License

This project is licensed under the Creative Commons Attribution 4.0
International License (CC BY 4.0) — see the [LICENSE](LICENSE) file for
details.

Please also cite the associated publication and the original data/database
sources listed under [References](#references) when using this workflow.
