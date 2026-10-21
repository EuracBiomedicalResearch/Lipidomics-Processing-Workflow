# Generate annotation metrics

Run from the repository root, using the project's Docker image or configured
R environment. Metrics do not change annotations, normalization or abundances.

## Refresh existing results

```bash
Rscript applications/MICROSAMPLING_study/generate_annotation_metrics.R
Rscript applications/METFORMIN-HIIE_study/generate_annotation_metrics.R
```

With the existing Docker image (run either study independently):

```bash
docker run --rm --user "$(id -u):$(id -g)" --mount "type=bind,source=$(pwd),target=/workspace" cembio-eurac:local MICROSAMPLING metrics
docker run --rm --user "$(id -u):$(id -g)" --mount "type=bind,source=$(pwd),target=/workspace" cembio-eurac:local METFORMIN metrics
```

Outputs are `applications/<STUDY>_study/objects/<STUDY>_annotation_metrics.xlsx`,
where `<STUDY>` is `MICROSAMPLING` or `METFORMIN-HIIE`.

Required inputs within each study:

- `positive/objects/<STUDY>_annotation_report_pos.rds`
- `negative/objects/<STUDY>_annotation_report_neg.rds`
- `objects/<STUDY>_annotation_report_merged.rds`

The first refresh also needs each polarity's saved
`<STUDY>_preprocessed_res_pos` / `<STUDY>_preprocessed_res_neg`, the study's
`LipidDatabase_R.xlsx`, and the corresponding `Annotation_pos.qmd` /
`Annotation_neg.qmd` to obtain the actual matching tolerances. It caches missing
mass-only and unadjusted-RT comparisons in the report RDS files. It does not
read raw spectra or rerun isotope validation. Later refreshes reuse these comparisons.

If reports are missing, run positive and negative annotation, then
`POS_NEG_merge.qmd`, reusing completed preprocessing. Merge automatically exports
metrics too. Missing/incompatible reports or inconsistent merge inputs cause an
error rather than reconstruction of final annotations. After changing matching
inputs or tolerances, rerun annotation and merge to create fresh snapshots;
cached comparisons are not automatically refreshed against changed inputs.

## Read the workbook

- `Summary`: actual workflow counts, positive then negative then merged.
- `Final_reference_comparison`: pair counts and three percentages against the
  fixed post-QC, pre-merge workflow reference, separately for each polarity.
- Remaining sheets: candidate/annotation rows for the recorded phases.

All sheets retain short descriptions below their column headers. ISTD labels
are excluded at every phase; assigned ISTD feature IDs (including standards
that fail QC) and their alternative labels are excluded throughout, even from
preprocessing counts. A mass-only standard candidate does not by itself identify
an entire feature as an ISTD. No internal-standards column is exported.

`ambiguity_resolution` is the display name for the legacy `manual_curation`
snapshot: automatic ambiguity decisions are applied and alternative labels are
consolidated. Its name does not assert human review. The underlying workflow
also permits reviewed Excel flags; an unattended render regenerates its default
ambiguity files before applying them.

For reference pairs G and candidates P_s, T_s counts pairs in both sets, C_s
counts alternative labels on reference features, and U_s counts pairs on other
features. They partition P_s. U0 is U_s at `rank1_mz` (rank-1 mass-only matching).

- Reference retention: `100 * T_s / |G|`.
- Non-reference candidate reduction: `100 * (1 - U_s / U0)`.
- Within-reference candidate discrimination: `100 * T_s / (T_s + C_s)`.

Pairs are distinct, polarity-qualified feature–lipid assignments; semicolon
alternatives count separately. Zero denominators display `N/A`. These are
workflow-reference comparisons, not independent chemical accuracy. Retention
along the filtering pipeline is normally 100% by construction.

`rank1_mz_rt_unadjusted` uses original database RTs with the same mass and RT
tolerances as adjusted matching. It is a parallel comparison, not a filtering
stage, so it appears in the comparison/detail sheets but not Summary.
