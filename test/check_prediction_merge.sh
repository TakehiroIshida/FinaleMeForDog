#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

PRED_DIR="out/pred_parts"
GLOB="${PRED_DIR}/prediction.part*.gz"
OUT_ALL="out/prediction.all.parts128.gz"

PRED_FILES=( $GLOB )
if [ "${#PRED_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no prediction files found: $GLOB" >&2
  exit 1
fi
if [ ! -f "$OUT_ALL" ]; then
  echo "ERROR: missing merged file: $OUT_ALL" >&2
  exit 1
fi

# sum of part lines
sum_lines=0
for f in "${PRED_FILES[@]}"; do
  l=$(zcat "$f" | wc -l)
  sum_lines=$((sum_lines + l))
done

merged_lines=$(zcat "$OUT_ALL" | wc -l)

# parts はそれぞれヘッダ1行ずつ持つので、結合時は (parts-1) 行減るのが正しい
expected=$((sum_lines - (${#PRED_FILES[@]} - 1)))

echo "INFO: parts files = ${#PRED_FILES[@]}" >&2
echo "INFO: sum(part lines) = $sum_lines" >&2
echo "INFO: merged lines    = $merged_lines" >&2
echo "INFO: expected merged = $expected" >&2

if [ "$merged_lines" -ne "$expected" ]; then
  echo "ERROR: merged line count mismatch" >&2
  exit 1
fi

echo "OK: merged line count matches expected" >&2

echo "INFO: merged head:" >&2
zcat "$OUT_ALL" | head -n 3 >&2
