#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

DETAILS="out/CpgMultiMetricsStats.G1000.low_cov.details.bed.gz"

OUT_DIR="out"
mkdir -p "$OUT_DIR"

# Usage:
#   ./make_input_matrix_full.sh [MAX_PER_READ]
#   MAX_PER_READ: 1 readName あたりの行数上限（OOM防止のため）。デフォルト20。
#
# Examples:
#   ./make_input_matrix_full.sh 20
#   ./make_input_matrix_full.sh 10
MAX_PER_READ="${1:-20}"

INPUT_MAT_FULL="${OUT_DIR}/input_matrix.details.full.max${MAX_PER_READ}.G1000.low_cov.tsv.gz"
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

# DETAILS (12列) -> INPUT_MAT_FULL (tsv.gz) へ変換
#  - FinaleMe.processMatrixFile と同じ行レベルのフィルタを適用:
#      fragLen: 30 < fragLen < 500
#      baseQ  : > 5
#      offset : >= 0
#      prior  : NaN/empty -> 0
#  - OOM（メモリ不足）回避のため readName ごとの出力行数を制限（MAX_PER_READ）
zcat "$DETAILS" \
| awk -v MAXR="$MAX_PER_READ" 'BEGIN{FS=OFS="\t"}
  NR==1 { print; next }

  # 12列ではない不要な行を除外
  NF!=12 { next }

  {
    fraglen = $5 + 0
    baseQ   = $9 + 0
    offset  = $10 + 0
    prior   = $12

    # FinaleMeと同じフィルタを適用
    if (fraglen >= 500 || fraglen <= 30) next
    if (baseQ <= 5) next
    if (offset < 0) next

    # prior が NaN/empty の場合は 0 に置換
    if (prior=="NaN" || prior=="nan" || prior=="") $12 = 0

    rn = $4

    # readName あたりの行数を制限（OOM防止）
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
