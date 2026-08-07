# CEMBIO-EURAC

**Automated data workflow for targeted lipid annotation using a public lipid database**

A vendor-independent, open-source workflow for LC-MS lipidomics data processing and annotation.

---

## Overview

This workflow provides:

- **Automated RT adjustment** using a Lipid reference set (containing SPLASH LIPIDOMIX and some endogenous compounds)
- **Targeted lipid annotation** against a curated SRM1950 human plasma lipid database
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

Each analysis is split by **ionization polarity** — positive and negative each
have their own Quarto documents so that parameter choices are fully traceable.

```
CEMBIO-EURAC/
│
├── generic_workflow/                # Reusable workflow templates
│   ├── LipidDatabase_R.xlsx         #   Lipid database (sheet 4=POS, sheet 5=NEG)
│   ├── positive/                    #   Positive-mode templates
│   │   ├── Preprocessing_pos.qmd    #     Step 1 (POLARITY = "pos")
│   │   └── Annotation_pos.qmd      #     Step 2 (POLARITY = "pos")
│   ├── negative/                    #   Negative-mode templates
│   │   ├── Preprocessing_neg.qmd    #     Step 1 (POLARITY = "neg")
│   │   └── Annotation_neg.qmd      #     Step 2 (POLARITY = "neg")
│   └── POS_NEG_merge.qmd           #   Step 3: Merge positive + negative results
│
├── applications/                    # Study-specific applications (self-contained)
│   ├── pilot_study/
│   │   ├── LipidDatabase_R.xlsx     #   Lipid database (shared by POS + NEG)
│   │   ├── positive/                #   Positive-mode analysis
│   │   │   ├── Preprocessing_pos.qmd   #     Configured for pilot, POS
│   │   │   ├── Annotation_pos.qmd      #     Configured for pilot, POS
│   │   │   ├── seq_pos_pilot.xlsx   #     Sample sequence
│   │   │   ├── pos_lipid_reference_set.xlsx  # Reference lipids
│   │   │   └── data/               #     .mzML files (gitignored)
│   │   ├── negative/                #   Negative-mode analysis
│   │   │   ├── Preprocessing_neg.qmd   #     Configured for pilot, NEG
│   │   │   ├── Annotation_neg.qmd      #     Configured for pilot, NEG
│   │   │   ├── seq_neg_pilot.xlsx   #     Sample sequence
│   │   │   ├── neg_lipid_reference_set.xlsx  # Reference lipids
│   │   │   └── data/               #     .mzML files (gitignored)
│   │   └── POS_NEG_merge.qmd       #   Step 3: Merge POS + NEG pilot results
│   └── exercise_study/              #   (same structure as pilot_study)
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

Each polarity folder is self-contained. Place these files inside your
polarity folder (e.g. `applications/my_study/positive/`):

1. **LC-MS raw data**: Place `.mzML` files in the `data/` subfolder

2. **Sample sequence file** (`seq_pos_<study_id>.xlsx` or `seq_neg_<study_id>.xlsx`):

   | file_name | sample_name | sample_type | injection_index |
   |-----------|-------------|-------------|------------------|
   | D01P_pos.mzML | D01P | Plasma | 1 |
   | QC_1_pos.mzML | QC_1 | QC | 2 |

3. **Lipid database**: `LipidDatabase_R.xlsx` in the study folder (one level above each polarity folder)

4. **Lipid Reference Set**: `pos_lipid_reference_set.xlsx` or `neg_lipid_reference_set.xlsx`

### 3. Run the Workflow

Each polarity has dedicated Quarto documents with hardcoded `POLARITY`,
so parameter choices are fully traceable.

| Step | Files | Description |
|------|-------|-------------|
| 1 | `positive/Preprocessing_pos.qmd`, `negative/Preprocessing_neg.qmd` | Data import, peak detection, RT alignment, gap filling |
| 2 | `positive/Annotation_pos.qmd`, `negative/Annotation_neg.qmd` | RT correction, database matching, normalization, QC |
| 3 | `POS_NEG_merge.qmd` | Combine POS + NEG results, coverage figures, PCA |

**Example: run the pilot study**

1. Render `applications/pilot_study/positive/Preprocessing_pos.qmd`
2. Render `applications/pilot_study/positive/Annotation_pos.qmd`
3. Render `applications/pilot_study/negative/Preprocessing_neg.qmd`
4. Render `applications/pilot_study/negative/Annotation_neg.qmd`
5. Render `applications/pilot_study/POS_NEG_merge.qmd`

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

## User Checkpoints

The workflow includes several interactive checkpoints:

| Checkpoint | Location | Action Required |
|------------|----------|-----------------|
| RT Filter Range | Preprocessing | Adjust RT filter based on BPC |
| Reference Lipid EICs | Preprocessing | Verify internal standard signals |
| Peak Detection Results | Preprocessing | Check peak detection quality |
| RT Correction Model | Annotation | Validate R-squared and residuals |
| Isotope Pattern Validation | Annotation | Review mirror plots |
| Manual Ambiguity Resolution | Annotation | Curate Excel files |

---

## Configuration Parameters

Key parameters in the Configuration section of each workflow file:

```r
PROJECT_ROOT <- "../../.."  # Path to project root (for shared R/ helper functions)
STUDY_ID <- "pilot"         # Study prefix for saved objects and seq file names
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

1. **Lipid Database**: <https://doi.org/10.1016/j.jlr.2024.100671>
2. **Original Study (pilot)**: <https://doi.org/10.1016/j.microc.2025.113760>
3. **Validation Study (exercise)**: <https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00480>

---

## Contributors

- CEMBIO-EURAC Team
- Sara Londono
- Philippe Louail

---

## License

This project is for research purposes. Please cite the original publications when using this workflow.
