#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

# ===== Java classpath =====
JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/*"

# ===== Inputs =====
MAX_PER_READ="${1:-20}"
INPUT_MAT="out/input_matrix.details.full.max${MAX_PER_READ}.tsv.gz"

# ===== Split settings =====
PARTS="${2:-128}"
OUT_DIR="out/parts"
LOG="out/split_input_matrix.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

# ===== Sanity checks =====
for f in "$JAR" "$INPUT_MAT"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

echo "INFO: INPUT_MAT=$INPUT_MAT" >&2
echo "INFO: PARTS=$PARTS" >&2
echo "INFO: OUT_DIR=$OUT_DIR" >&2

# ===== Run splitter =====
# Expected signature:
#   SplitInputMatrix <input.tsv.gz> <out_dir> <parts>
java -Xmx4G -cp "$CP" \
  org.cchmc.epifluidlab.finaleme.utils.SplitInputMatrix \
  "$INPUT_MAT" "$OUT_DIR" "$PARTS" \
  2> "$LOG"

echo "OK: split completed" >&2
echo "LOG: $LOG" >&2

# ===== Quick checks =====
echo "INFO: counting lines..." >&2

TOTAL_ORIG=$(zcat "$INPUT_MAT" | wc -l)

TOTAL_PARTS=0
PART_FILES=( "$OUT_DIR"/input_matrix.part*.tsv.gz )
if [ "${#PART_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no part files created in $OUT_DIR" >&2
  exit 1
fi

for f in "${PART_FILES[@]}"; do
  LINES=$(zcat "$f" | wc -l)
  TOTAL_PARTS=$((TOTAL_PARTS + LINES))
done

echo "INFO: original lines = $TOTAL_ORIG" >&2
echo "INFO: parts sum lines = $TOTAL_PARTS" >&2

if [ "$TOTAL_ORIG" -ne "$TOTAL_PARTS" ]; then
  echo "ERROR: line count mismatch (orig=$TOTAL_ORIG, parts=$TOTAL_PARTS)" >&2
  exit 1
fi

echo "OK: line counts match" >&2

# Optional: show top 5 largest parts by line count
echo "INFO: top 5 largest parts (by lines):" >&2
for f in "${PART_FILES[@]}"; do
  echo -e "$(zcat "$f" | wc -l)\t$f"
done | sort -nr | head -n 5 >&2
