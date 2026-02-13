#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

PRED_DIR="out/pred_parts"
GLOB="${PRED_DIR}/prediction.part*.gz"

OUT_DIR="out"
mkdir -p "$OUT_DIR"

OUT_ALL="${OUT_DIR}/prediction.all.gz"
TMP_OUT="${OUT_ALL}.tmp.$$"
rm -f "$TMP_OUT"

PRED_FILES=( $GLOB )
if [ "${#PRED_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no prediction files found: $GLOB" >&2
  exit 1
fi

IFS=$'\n' PRED_FILES_SORTED=($(printf "%s\n" "${PRED_FILES[@]}" | sort))
unset IFS

echo "INFO: merging ${#PRED_FILES_SORTED[@]} files" >&2

# ヘッダは python で安全に取得
HEADER="$(python - << 'PY'
import gzip, sys
import glob
files = sorted(glob.glob("out/pred_parts/prediction.part*.gz"))
with gzip.open(files[0], "rt") as f:
    print(f.readline().rstrip("\n"))
PY
)"

{
  printf '%s\n' "$HEADER"
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
