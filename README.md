# CEMBIO-EURAC

**Automated data workflow for targeted lipid annotation using a public lipid database**

A vendor-independent, open-source workflow for LC-MS lipidomics data processing and annotation.

---

## 📋 Overview

This workflow provides:

- **Automated RT correction** using reference lipids (SPLASH LIPIDOMIX)
- **Targeted lipid annotation** against a curated human plasma lipid database
- **Isotope pattern validation** for high-confidence identification
- **Adduct profile matching** for improved annotation
- **Internal standard normalization** by lipid subclass
- **Quality control filtering** based on QC sample RSD

### Supported Data

- **Ionization modes**: Positive (POS) and Negative (NEG)
- **Instrument platforms**: Vendor-independent (requires .mzML format)
- **Sample types**: QC, M-VAMS, W-DBS, C-qDBS, Plasma (configurable)

---

## 📁 Project Structure

```
CEMBIO-EURAC/
│
├── Lipidomics_workflow.qmd      # 📌 MAIN WORKFLOW - Start here!
│
├── R/
│   └── lipid_helpers.R          # Helper functions (loaded automatically)
│
├── POS_data/                    # Positive mode raw data (.mzML files)
│   └── *.mzML
│
├── NEG_data/                    # Negative mode raw data (.mzML files)
│   └── *.mzML
│
├── objects/                     # Saved R objects (auto-generated)
│   ├── preprocessed_mse_pos/    # Preprocessed MsExperiment (positive)
│   ├── preprocessed_mse_neg/    # Preprocessed MsExperiment (negative)
│   ├── preprocessed_res/        # SummarizedExperiment (positive)
│   └── preprocessed_res_neg/    # SummarizedExperiment (negative)
│
├── figures/                     # Output figures (auto-generated)
│   ├── EIC_internal_standards/  # EIC plots for reference lipids
│   ├── iso_pattern_check/       # Isotope pattern mirror plots (pos)
│   ├── iso_pattern_check_neg/   # Isotope pattern mirror plots (neg)
│   ├── ref_lipid_image/         # RT correction diagnostics (pos)
│   └── ref_lipid_image_neg/     # RT correction diagnostics (neg)
│
├── pos_peak_detection_ref_lipid/  # Peak detection QC plots (pos)
├── neg_peak_detection_ref_lipid/  # Peak detection QC plots (neg)
├── ISTD_mtched_data/              # Matched ISTD chromatograms (pos)
├── ISTD_mtched_data_neg/          # Matched ISTD chromatograms (neg)
│
├── seq_pos.xlsx                 # Sample sequence (positive mode)
├── seq_neg.xlsx                 # Sample sequence (negative mode)
├── splashlipidomix_list.xlsx    # Reference lipids (positive mode)
├── neg_splashlipidomix_list.xlsx # Reference lipids (negative mode)
├── LipidDatabase_R.xlsx         # Lipid database (sheet 4=POS, sheet 5=NEG)
│
├── Pilot_POS.qmd                # Original positive mode analysis
├── Pilot_NEG.qmd                # Original negative mode analysis
└── README.md                    # This file
```

---

## 🚀 Quick Start

### 1. Prerequisites

Install required R packages:

```r
# Bioconductor packages
BiocManager::install(c(
  "MsExperiment", "MsIO", "alabaster.se", "MsBackendMetaboLights",
  "SummarizedExperiment", "xcms", "Spectra", "MetaboCoreUtils",
  "limma", "matrixStats", "BiocFileCache", "AnnotationHub",
  "CompoundDb", "MetaboAnnotation"
))

# CRAN packages
install.packages(c(
  "knitr", "readxl", "writexl", "pander", "RColorBrewer",
  "pheatmap", "vioplot", "ggplot2", "ggfortify", "gridExtra",
  "enviPat", "ggVennDiagram", "UpSetR", "dbplyr"
))
```

### 2. Prepare Input Files

1. **Sample sequence file** (`seq_pos.xlsx` or `seq_neg.xlsx`):
   | file_name | sample_name | sample_type | injection_index |
   |-----------|-------------|-------------|-----------------|
   | D01P_pos.mzML | D01P | Plasma | 1 |
   | QC_1_pos.mzML | QC_1 | QC | 2 |
   | ... | ... | ... | ... |

2. **Reference lipid list**: Ensure your internal standards file has columns `short_name`, `mz`, `RT` (in seconds)

3. **Raw data**: Convert your vendor files to `.mzML` format and place in `POS_data/` or `NEG_data/`

### 3. Run the Workflow

1. Open `Lipidomics_workflow.qmd` in RStudio/Positron
2. Set the polarity in the Configuration section:
   ```r
   POLARITY <- "pos"  # or "neg"
   ```
3. Adjust parameters as needed (RT range, peak detection, etc.)
4. Render the document or run chunks interactively

---

## 📋 User Checkpoints

The workflow includes several **USER CHECKPOINT** sections that require your attention:

| Checkpoint | Location | Action Required |
|------------|----------|-----------------|
| RT Filter Range | Data Import | Adjust RT filter based on BPC |
| Reference Lipid EICs | Dataset Investigation | Verify IS signals |
| Peak Detection Results | Preprocessing | Check peak detection quality |
| RT Correction Model | Database Correction | Validate R² and residuals |
| Isotope Pattern Validation | Annotation | Review mirror plots |
| Manual Ambiguity Resolution | Annotation | Curate Excel files |

---

## 📊 Output Files

### Generated During Analysis

| File | Description |
|------|-------------|
| `lipid_ambiguity_resolution.xlsx` | Lipids matching multiple features (manual curation) |
| `feature_ambiguity_resolution.xlsx` | Features matching multiple lipids (manual curation) |

### Saved R Objects

| Object | Contents |
|--------|----------|
| `preprocessed_mse_*` | MsExperiment with chromatographic peaks |
| `preprocessed_res*` | SummarizedExperiment with feature intensities |

---

## 🔧 Configuration Parameters

Key parameters in `Lipidomics_workflow.qmd`:

```r
POLARITY <- "pos"          # Analysis mode
CORES_NB <- 4              # Parallel processing cores
RT_FILTER_MIN <- 10        # RT filter start (seconds)
RT_FILTER_MAX <- 950       # RT filter end (seconds)
PEAK_WIDTH <- c(4, 8)      # Peak width range (seconds)
PPM <- 10                  # m/z tolerance for peak detection
MATCH_PPM <- 20            # m/z tolerance for database matching
MATCH_RT_TOL <- 20         # RT tolerance for matching (seconds)
ISOPEAK_SIM_THRESHOLD <- 0.78  # Minimum isotope similarity
RSD_THRESHOLD <- 0.3       # QC RSD filter (30%)
```

---

## 📚 References

1. **Lipid Database**: [https://doi.org/10.1016/j.jlr.2024.100671](https://doi.org/10.1016/j.jlr.2024.100671)
2. **Original Study**: [https://doi.org/10.1016/j.microc.2025.113760](https://doi.org/10.1016/j.microc.2025.113760)

---

## 👥 Contributors

- CEMBIO-EURAC Team
- Sara Londono
- Philippe Louail

---

## 📄 License

This project is for research purposes. Please cite the original publications when using this workflow.
