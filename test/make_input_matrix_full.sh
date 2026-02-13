#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

DETAILS="out/CpgMultiMetricsStats.hg19.details.bed.gz"

OUT_DIR="out"
mkdir -p "$OUT_DIR"

# Usage:
#   ./make_input_matrix_full.sh [MAX_PER_READ]
#
# Examples:
#   ./make_input_matrix_full.sh 20
#   ./make_input_matrix_full.sh 10
MAX_PER_READ="${1:-20}"

INPUT_MAT_FULL="${OUT_DIR}/input_matrix.details.full.max${MAX_PER_READ}.tsv.gz"
TMP_OUT="${INPUT_MAT_FULL}.tmp.$$"
rm -f "$TMP_OUT"

for f in "$DETAILS"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

echo "INFO: DETAILS=$DETAILS" >&2
echo "INFO: MAX_PER_READ=$MAX_PER_READ" >&2
echo "INFO: OUTPUT=$INPUT_MAT_FULL" >&2

# DETAILS (12 cols) -> INPUT_MAT_FULL (tsv.gz)
#  - Keep the same row-level filters as FinaleMe.processMatrixFile:
#      fragLen: 30 < fragLen < 500
#      baseQ  : > 5
#      offset : >= 0
#      prior  : NaN/empty -> 0
#  - Additionally limit number of lines per readName (MAX_PER_READ) to avoid OOM.
zcat "$DETAILS" \
| awk -v MAXR="$MAX_PER_READ" 'BEGIN{FS=OFS="\t"}
  NR==1 { print; next }

  # Drop malformed lines
  NF!=12 { next }

  {
    fraglen = $5 + 0
    baseQ   = $9 + 0
    offset  = $10 + 0
    prior   = $12

    # Align with FinaleMe row-level filters
    if (fraglen >= 500 || fraglen <= 30) next
    if (baseQ <= 5) next
    if (offset < 0) next

    # Prior cleanup
    if (prior=="NaN" || prior=="nan" || prior=="") $12 = 0

    rn = $4

    # Cap lines per fragment/readName (OOM guard)
    if (++seen[rn] > MAXR) next

    print
  }
  END{
    total=0
    for (x in seen) total += seen[x]
    print "INFO: total lines written =", total > "/dev/stderr"
    print "INFO: total fragments observed =", length(seen) > "/dev/stderr"
  }' \
| gzip -c > "$TMP_OUT"

gzip -t "$TMP_OUT"
mv -f "$TMP_OUT" "$INPUT_MAT_FULL"

echo "OK: wrote $INPUT_MAT_FULL" >&2
ls -lh "$INPUT_MAT_FULL" >&2
echo -n "INFO: INPUT_MAT_FULL lines=" >&2
zcat "$INPUT_MAT_FULL" | wc -l >&2
