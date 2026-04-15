# PosNegExp Extension for Annotation — Design Analysis

## Context

`applications/pilot_study/Annotation.qmd` currently uses explicit
`for (POLARITY in c("pos", "neg"))` loops to apply each annotation,
normalization, and QC step to both ionization modes. This document analyses
which loops could benefit from a `PosNegExp`-style dispatch approach — and
where a simpler or different strategy is more appropriate.

---

## Why PosNegExp dispatch works well for Preprocessing

In `Preprocessing.qmd`, the dispatched functions are **xcms S4 generics**
(`findChromPeaks`, `adjustRtime`, `groupChromPeaks`, `fillChromPeaks`).
Dispatch on a new class is natural here because:

1. xcms already defines these as `setGeneric()` — adding a `setMethod()` for
   `PosNegXcmsExp` follows the established pattern.
2. `PosNegXcmsExp` fits cleanly into the Bioconductor class hierarchy.
3. The user-facing API is unchanged — they call the same xcms functions.
4. Parallelisation and memory handling are hidden inside the dispatch.

---

## Why the same approach is also right for Annotation

Two arguments justify extending the wrapper-class pattern into annotation,
beyond any per-chunk cost-benefit calculation:

**1. Custom functions do not need full S4 dispatch.**
The original analysis assumed that making a function aware of a wrapper class
requires `setGeneric()` + `setMethod()` boilerplate (8–10 lines per function).
That is true for xcms generics we do not own. For our own functions in
`lipid_helpers.R` the cost is an `inherits()` guard at the function head:

```r
apply_volume_correction <- function(se, factors, assay_name, new_assay_name) {
    if (inherits(se, "PosNegSumExp"))
        return(PosNegSumExp(
            apply_volume_correction(se@pos, factors, assay_name, new_assay_name),
            apply_volume_correction(se@neg, factors, assay_name, new_assay_name)
        ))
    # existing single-mode implementation unchanged...
}
```

That is ~4 lines per function. Three Tier 1 functions → 12 lines of guards,
plus ~20 lines for the class definition. No new generics, no interface
contracts.

**2. Consistency with Preprocessing is a design requirement.**
Preprocessing establishes the paradigm: *wrap → call dispatched function →
result stays wrapped*. Annotation immediately follows Preprocessing in the
workflow. If annotation breaks back to explicit `for (POLARITY ...)` loops
throughout, readers face a paradigm shift in the same pipeline. Where the
wrapper pattern fits, it should be used.

---

## Annotation functions: three tiers

### Tier 1 — Pure `SummarizedExperiment` transformers (identical call, both modes)

| Chunk | Function | Signature |
|---|---|---|
| `volume-correction` | `apply_volume_correction` | `(se, factors, assay_name, new_assay_name)` |
| `is-normalization` | `normalize_by_is` | `(se, input_assay, output_assay, is_col, lipid_col)` |
| `rsd-filter` | `filter_by_qc_rsd` | `(se, threshold, qc_col, qc_value)` |

These three take a `SummarizedExperiment` in and return one out with identical
args for both modes. Since all three are custom functions in `R/lipid_helpers.R`,
the cost is an `inherits()` guard (not full S4 dispatch — see above).

**`PosNegSumExp` class** (~20 lines in `R/pos_neg_exp.R`):

```r
setClass("PosNegSumExp",
    slots = c(pos = "SummarizedExperiment", neg = "SummarizedExperiment"))
PosNegSumExp <- function(pos, neg) new("PosNegSumExp", pos = pos, neg = neg)
posRes <- function(x) x@pos
negRes <- function(x) x@neg
```

With this, the `volume-correction` loop collapses from:

```r
for (POLARITY in c("pos", "neg")) {
    res <- if (POLARITY == "pos") res_pos else res_neg
    res <- apply_volume_correction(res, VOLUME_FACTORS,
                                   assay_name = "raw", new_assay_name = "raw_corr")
    res <- apply_volume_correction(res, VOLUME_FACTORS,
                                   assay_name = "raw_filled", new_assay_name = "corr_filled")
    res_corr_list[[POLARITY]] <- res
}
```

to:

```r
res <- apply_volume_correction(res, VOLUME_FACTORS,
                               assay_name = "raw", new_assay_name = "raw_corr")
res <- apply_volume_correction(res, VOLUME_FACTORS,
                               assay_name = "raw_filled", new_assay_name = "corr_filled")
```

The same pattern applies to `normalize_by_is` and `filter_by_qc_rsd`.

**Plotting loops (`pca-plot`, `rla-plot`, `cv-summary`) stay as explicit loops.**
They require per-mode titles, per-mode color data, and print separate figures.
Even with `PosNegSumExp` you would unwrap inside the loop — no gain.

---

### Tier 2 — Functions with mode-specific arguments

| Chunk | Function | Differing args |
|---|---|---|
| `validate-database` | `validate_lipid_database` | `sheet`, `polarity`, `rt_col` |
| `rt-correction-fit` | `fit_rt_correction` | `eic_is`, `intern_standard` |
| `load-database` | `prepare_lipid_database` | `sheet`, `polarity`, `rt_col`, `rt_fit` |
| `plot-rt-correction` | inline plot | `rt_fit`, `lipid_database`, `intern_standard` |
| `rank1-matching` | `match_features_to_database` | `res`, `lipid_database` |
| `isotope-validation` | `calculate_isotope_similarity` | `mse`, `polarity` |

These all have genuinely per-mode inputs or outputs. The RT correction fit
produces a mode-specific polynomial model; database loading yields a different
database per mode; rank 1 matching consumes both — making them candidates for
loop retention.

**`calculate_isotope_similarity`** is the exception: it bridges the xcms world
(per-mode `mse_pos`/`mse_neg`) and the annotation world. A single new
`setMethod("calculate_isotope_similarity", "PosNegXcmsExp", ...)` that
internally extracts `posExp()`/`negExp()` and dispatches per mode is
non-breaking and fits the existing `PosNegExp` pattern. Keep the plot loop
separate.

**All other Tier 2 loops stay explicit** — their outputs are necessarily
per-mode and there is no meaningful type to dispatch on that would avoid
passing mode-specific arguments.

---

### Tier 3 — Annotation state: `AnnotationResult`

| Chunk | Core function | Primary argument |
|---|---|---|
| `adduct-matching` | `match_adducts` | `mtched_data` + `lipid_database` + `query` |
| `resolve-ambiguities` | inline deduplication logic | `mtched_data` |
| `sm-resolution` | `resolve_sm_isomers` | `mtched_data` |
| `export-ambiguity` | `export_ambiguity_tables` | `mtched_data` + per-mode filenames |
| `annotation-summary` | `plot_lipid_donut` | `mtched_data` |

Every function in this tier operates on the same `mtched_data` schema, and
`match_adducts` always needs three pieces of context together:

```r
match_adducts(mtched_data, lipid_database, query, ppm, rt_tol)
```

These three always travel as a group from `match_features_to_database` onward.
The classic sign that a class is warranted.

#### The right container: `AnnotationResult`

```r
setClass("AnnotationResult",
    slots = c(
        matches  = "DataFrame",   # mtched_data — the evolving annotation state
        query    = "DataFrame",   # experimental features from res
        database = "DataFrame"    # full lipid_database for this mode
    )
)
```

Since these are our own functions, dispatch uses `inherits()` guards (same
rationale as Tier 1):

```r
match_adducts <- function(mtched_data, lipid_database, query, ppm, rt_tol) {
    if (inherits(mtched_data, "AnnotationResult"))
        return(AnnotationResult(
            matches  = match_adducts(mtched_data@matches, mtched_data@database,
                                     mtched_data@query, ppm, rt_tol),
            query    = mtched_data@query,
            database = mtched_data@database
        ))
    # existing implementation...
}
```

`resolve_sm_isomers` and `export_ambiguity_tables` reduce to single-argument
dispatch the same way.

#### Prerequisite: extract `resolve-ambiguities` inline logic

The `resolve-ambiguities` chunk contains ~30 lines of key-counting and
RT-score-filtering logic written directly in the qmd. This must become a
function `resolve_annotation_ambiguities(annotation)` in `lipid_helpers.R`
before it can participate in dispatch. This extraction is required regardless
of whether `AnnotationResult` is adopted.

#### `PosNegAnnotation` class

```r
setClass("PosNegAnnotation",
    slots = c(pos = "AnnotationResult", neg = "AnnotationResult"))
```

With `inherits()` guards extended to handle `PosNegAnnotation`, the
adduct-matching through summary chain collapses to:

```r
annotation <- match_adducts(annotation, ppm = MATCH_PPM, rt_tol = 5)
annotation <- resolve_annotation_ambiguities(annotation)
annotation <- resolve_sm_isomers(annotation)
export_ambiguity_tables(annotation,
    lipid_file   = perMode(pos = "lipid_ambiguity_resolution_pos.xlsx",
                           neg = "lipid_ambiguity_resolution_neg.xlsx"),
    feature_file = perMode(pos = "feature_ambiguity_resolution_pos.xlsx",
                           neg = "feature_ambiguity_resolution_neg.xlsx"))
```

The `apply-curation` loop stays explicit — it reads per-mode xlsx files,
applies manual curation decisions, and merges duplicate features. This logic
is inherently per-mode and re-runnability depends on it being transparent.

---

### The missing bridge: `annotate_features()`

The `annotated-res` chunk joins the SE world (`res_corr_list`) with the
annotation world (`mtched_data_list`):

```r
annotated_res <- res[rownames(res) %in% mtched_data$feature_id, ]
rowData(annotated_res) <- cbind(rowData(annotated_res), mtched_data[...])
```

With both `PosNegSumExp` and `PosNegAnnotation` in place, this step needs a
**bridge function** that is not currently in the helpers:

```r
annotate_features <- function(res, annotation) {
    # res: PosNegSumExp; annotation: PosNegAnnotation
    # returns: PosNegSumExp with rowData populated from annotation
    PosNegSumExp(
        .annotate_se(posRes(res), annotation@pos),
        .annotate_se(negRes(res), annotation@neg)
    )
}
```

Without this function, the `annotated-res` construction loop cannot collapse
even if both wrapper types exist.

---

## Recommendation

### Three extension points, ordered by value

#### 1. `AnnotationResult` + `PosNegAnnotation` — highest value

The `mtched_data` + `query` + `database` triplet is a de-facto object with a
stable schema. Formalising it eliminates the longest and most complex loops in
the qmd. The `inherits()` guard approach keeps the cost proportional.

#### 2. `PosNegSumExp` — moderate value, low cost

Adds consistency with Preprocessing for the SE-transform loops. Cost is ~30
lines total (class + guards in 3 functions). Not a paradigm shift — just an
`inherits()` guard in each custom function.

#### 3. `calculate_isotope_similarity` `PosNegXcmsExp` method — contained, fits pattern

Single new `setMethod`, no signature changes, non-breaking.

#### Keep explicit loops for:

- RT correction fit, database loading, RT diagnostic plot — per-mode outputs
- Validation and data loading — per-mode file paths
- `apply-curation` — per-mode xlsx paths, re-runnability logic
- Plotting loops (`pca-plot`, `rla-plot`, `cv-summary`, `remove-qc`) — per-mode
  visualization with distinct titles, colors, and file paths

---

## Summary table

| Chunk group | Current | Recommended |
|---|---|---|
| Validate DB, load preprocessed | `for` loop | keep loop (per-mode paths/args) |
| RT correction fit, DB load, RT plot | `for` loop | keep loop (per-mode outputs) |
| Rank 1 matching | `for` loop | keep loop; update to return `AnnotationResult` |
| Isotope validation | `for` loop | `PosNegXcmsExp` method on `calculate_isotope_similarity`; keep plot loop |
| Adduct matching, ambiguity resolution, SM | `for` loop | **`AnnotationResult` + `PosNegAnnotation` dispatch** |
| Export ambiguity | `for` loop | `PosNegAnnotation` method with `perMode()` filenames |
| Apply curation | `for` loop | keep loop (per-mode xlsx, re-runnability) |
| Annotation summary | `for` loop | **`PosNegAnnotation` dispatch** |
| Volume correction | `for` loop | **`PosNegSumExp` + `inherits()` guard** |
| Annotated-res construction | `for` loop | **`annotate_features()` bridge** |
| Imputation, IS norm, RSD filter | `for` loop | **`PosNegSumExp` + `inherits()` guards** |
| PCA, RLA, CV plots, remove-QC, save | `for` loop | keep loop (per-mode visualization/I/O) |

## Implementation order

1. **Extract `resolve-ambiguities` inline logic** → `resolve_annotation_ambiguities()`
   in `R/lipid_helpers.R`. Required before any dispatch on the annotation chain.

2. **`AnnotationResult` class** in `R/pos_neg_exp.R` or a new
   `R/annotation_result.R` — define schema, constructor, accessors
   (`matches()`, `query()`, `database()`). Update `match_features_to_database`
   to return one.

3. **Add `inherits()` guards** to `match_adducts`, `resolve_annotation_ambiguities`,
   `resolve_sm_isomers`, `export_ambiguity_tables`, `plot_lipid_donut`
   dispatching on `AnnotationResult`.

4. **`PosNegAnnotation` class** — wrap two `AnnotationResult` objects; extend
   the same guards to handle `PosNegAnnotation` and collapse the annotation
   chain loops.

5. **`PosNegSumExp` class** (~20 lines in `R/pos_neg_exp.R`) + `inherits()`
   guards in `apply_volume_correction`, `normalize_by_is`, `filter_by_qc_rsd`.

6. **`annotate_features()` bridge function** — joins `PosNegSumExp` ×
   `PosNegAnnotation` → `PosNegSumExp`. Required for the `annotated-res`
   construction step to collapse.

7. **`calculate_isotope_similarity` `PosNegXcmsExp` method** — single new
   `setMethod`, no signature changes.
