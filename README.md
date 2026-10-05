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

### 1. Prerequisites

Install required R packages:

```r
# Bioconductor packages
BiocManager::install(c(
  "MsExperiment", "alabaster.se", "MsBackendMetaboLights",
  "SummarizedExperiment", "xcms", "Spectra", "MetaboCoreUtils",
  "limma", "matrixStats", "BiocFileCache", "AnnotationHub",
  "CompoundDb", "MetaboAnnotation"
))

# MsIO — pin to version 0.0.15
remotes::install_version("MsIO", version = "0.0.15")

# CRAN packages
install.packages(c(
  "knitr", "readxl", "writexl", "pander", "RColorBrewer",
  "pheatmap", "vioplot", "ggplot2", "ggfortify", "gridExtra",
  "enviPat", "ggVennDiagram", "UpSetR", "dbplyr"
))
```

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

**Example: run the MICROSAMPLING study**

1. Render `applications/MICROSAMPLING_study/positive/Preprocessing_pos.qmd`
2. Render `applications/MICROSAMPLING_study/positive/Annotation_pos.qmd`
3. Render `applications/MICROSAMPLING_study/negative/Preprocessing_neg.qmd`
4. Render `applications/MICROSAMPLING_study/negative/Annotation_neg.qmd`
5. Render `applications/MICROSAMPLING_study/POS_NEG_merge.qmd`

#### Annotation metrics

Both studies capture annotation-stage snapshots and export
`objects/<STUDY_ID>_annotation_metrics.xlsx` at the end of the merge.
Install the additional reporting dependency with `install.packages("openxlsx")`.
Existing preprocessing outputs can be reused, but rerun both annotation documents
and the merge once to capture the new snapshots. Missing or incompatible snapshots
produce an error instead of reconstructing annotations with separate settings.

The workbook starts with `Summary`: actual workflow phases, positive then negative
then merged, with feature counts, distinct feature–lipid pairs, removed features,
and ambiguous features. Standards remain in intermediate counts and details until
the explicit final removal phase. The final detail sheet matches the feature rows
in the annotated-abundances workbook. Column definitions appear beneath headers.
`Curated_reference_comparison` reports agreement with manually curated assignments,
not accuracy against independently established chemical identities.

After snapshots exist, refresh reporting without rerunning annotation:

```bash
Rscript applications/MICROSAMPLING_study/generate_annotation_metrics.R
Rscript applications/METFORMIN-HIIE_study/generate_annotation_metrics.R
```

Reporting regression tests: `Rscript tests/test_annotation_reporting.R`.

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
