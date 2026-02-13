#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

PRED_DIR="out/pred_parts"
GLOB="${PRED_DIR}/prediction.part*.gz"

OUT_DIR="out"
mkdir -p "$OUT_DIR"

OUT_ALL="${OUT_DIR}/prediction.all.parts128.gz"
TMP_OUT="${OUT_ALL}.tmp.$$"
rm -f "$TMP_OUT"

PRED_FILES=( $GLOB )
if [ "${#PRED_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no prediction files found: $GLOB" >&2
  exit 1
fi

# sort by part number (lexicographic works with 000..127)
IFS=$'\n' PRED_FILES_SORTED=($(printf "%s\n" "${PRED_FILES[@]}" | sort))
unset IFS

echo "INFO: merging ${#PRED_FILES_SORTED[@]} files into $OUT_ALL" >&2

# header は最初のファイルだけ採用し、以降は2行目以降を連結
{
  zcat "${PRED_FILES_SORTED[0]}" | head -n 1
  for f in "${PRED_FILES_SORTED[@]}"; do
    zcat "$f" | awk 'NR>1'
  done
} | gzip -c > "$TMP_OUT"

gzip -t "$TMP_OUT"
mv -f "$TMP_OUT" "$OUT_ALL"

echo "OK: wrote $OUT_ALL" >&2
ls -lh "$OUT_ALL" >&2
echo -n "INFO: OUT_ALL lines=" >&2
zcat "$OUT_ALL" | wc -l >&2
