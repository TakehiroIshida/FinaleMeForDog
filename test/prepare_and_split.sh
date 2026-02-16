#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

DETAILS="out/CpgMultiMetricsStats.hg19.details.bed.gz"

PARTS="${1:-128}"
MAX_PER_READ="${2:-20}"
OUT_DIR="out/parts"
LOG="out/prepare_and_split.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="target/classes:${JAR}:lib/*"

for f in "$DETAILS"; do
  [ -f "$f" ] || { echo "ERROR: missing file: $f" >&2; exit 1; }
done

echo "INFO: DETAILS=$DETAILS" >&2
echo "INFO: PARTS=$PARTS" >&2
echo "INFO: MAX_PER_READ=$MAX_PER_READ" >&2
echo "INFO: OUT_DIR=$OUT_DIR" >&2
echo "INFO: LOG=$LOG" >&2

mvn -q -DskipTests compile

java -Xmx6G -cp "$CP" \
  org.cchmc.epifluidlab.finaleme.utils.PrepareAndSplitInputMatrix \
  "$DETAILS" "$OUT_DIR" "$PARTS" "$MAX_PER_READ" \
  2> "$LOG"

echo "OK: prepare+split completed" >&2

# ファイル数を確認
ls -1 "$OUT_DIR"/input_matrix.part*.tsv.gz | wc -l >&2

# 各partの行数を確認（上位5件）
for f in "$OUT_DIR"/input_matrix.part*.tsv.gz; do
  echo -e "$(zcat "$f" | wc -l)\t$f"
done | sort -nr | head -n 5 >&2
