#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

# ===== Inputs =====
MAX_PER_READ="${1:-20}"
INPUT_MAT="out/input_matrix.details.full.max${MAX_PER_READ}.tsv.gz"

# ===== Split settings =====
PARTS="${2:-128}"
OUT_DIR="out/parts"
LOG="out/split_input_matrix.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

# ===== Build/Classpath =====
# 重要: 新規作成したクラスは jar に入ってない可能性があるので target/classes を足す
# Mavenの場合、事前に mvn -q -DskipTests package か mvn -q -DskipTests compile を推奨
JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="target/classes:${JAR}:lib/*"

# ===== Sanity checks =====
for f in "$INPUT_MAT"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

if [ ! -f "$JAR" ]; then
  echo "WARN: missing $JAR (ok if you rely on target/classes + lib/* only)" >&2
fi

echo "INFO: INPUT_MAT=$INPUT_MAT" >&2
echo "INFO: PARTS=$PARTS" >&2
echo "INFO: OUT_DIR=$OUT_DIR" >&2
echo "INFO: LOG=$LOG" >&2

# ===== Run splitter =====
# 正しいFQCNに注意（tools）
java -Xmx4G -cp "$CP" \
  org.cchmc.epifluidlab.finaleme.tools.SplitInputMatrix \
  "$INPUT_MAT" "$OUT_DIR" "$PARTS" \
  2> "$LOG"

echo "OK: split completed" >&2

# ===== Quick checks =====
PART_FILES=( "$OUT_DIR"/input_matrix.part*.tsv.gz )
if [ "${#PART_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no part files created in $OUT_DIR" >&2
  echo "---- tail log ----" >&2
  tail -n 80 "$LOG" >&2 || true
  exit 1
fi

TOTAL_ORIG=$(zcat "$INPUT_MAT" | wc -l)
TOTAL_PARTS=0
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

echo "INFO: top 5 largest parts (by lines):" >&2
for f in "${PART_FILES[@]}"; do
  echo -e "$(zcat "$f" | wc -l)\t$f"
done | sort -nr | head -n 5 >&2
