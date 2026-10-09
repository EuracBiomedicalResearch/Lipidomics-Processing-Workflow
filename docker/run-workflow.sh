#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run-workflow STUDY [STAGE] [POLARITY]
  STUDY:     MICROSAMPLING | METFORMIN
  STAGE:     all (default) | preprocessing | annotation | merge
  POLARITY:  both (default) | pos | neg

Examples:
  run-workflow MICROSAMPLING
  run-workflow METFORMIN preprocessing pos
  run-workflow MICROSAMPLING annotation
  run-workflow MICROSAMPLING merge

Other commands: help | check | bash | Rscript ... | quarto ...
Run from the repository root; Quarto executes each document in its own folder.
EOF
}

case "${1:-help}" in
    help|-h|--help) usage; exit 0 ;;
    check) cd /workspace; exec Rscript scripts/check_environment.R ;;
    MICROSAMPLING|MICROSAMPLING_study) study=MICROSAMPLING_study ;;
    METFORMIN|METFORMIN-HIIE|METFORMIN-HIIE_study) study=METFORMIN-HIIE_study ;;
    *) exec "$@" ;;
esac

stage=${2:-all}
polarity=${3:-both}
if (( $# > 3 )); then usage >&2; exit 2; fi
case "$stage" in
    all|preprocessing|annotation|merge) ;;
    *) echo "Unknown stage: $stage" >&2; usage >&2; exit 2 ;;
esac
case "$polarity" in
    both) polarities=(pos neg) ;;
    pos|neg) polarities=("$polarity") ;;
    *) echo "Unknown polarity: $polarity" >&2; usage >&2; exit 2 ;;
esac
if [[ "$polarity" != both && "$stage" =~ ^(all|merge)$ ]]; then
    echo "$stage requires both polarities; omit POLARITY or use 'both'." >&2
    exit 2
fi

cd /workspace
mkdir -p "${BFC_CACHE:-/cache/BiocFileCache}"
study_dir="applications/$study"
render() {
    local document=$1
    echo "Rendering $document"
    # Explicit cwd preserves relative input paths and serialized backend paths.
    (cd "$(dirname "$document")"; quarto render "$(basename "$document")")
}

for mode in "${polarities[@]}"; do
    folder=positive
    [[ "$mode" == neg ]] && folder=negative
    if [[ "$stage" == all || "$stage" == preprocessing ]]; then
        render "$study_dir/$folder/Preprocessing_$mode.qmd"
    fi
    if [[ "$stage" == all || "$stage" == annotation ]]; then
        render "$study_dir/$folder/Annotation_$mode.qmd"
    fi
done
if [[ "$stage" == all || "$stage" == merge ]]; then
    render "$study_dir/POS_NEG_merge.qmd"
fi
