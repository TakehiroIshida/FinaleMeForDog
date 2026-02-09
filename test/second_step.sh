#!/usr/bin/env bash
set -euo pipefail

# --- move to project root (FinaleMeForDog/) even if run from test/ ---
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

# --- inputs / outputs ---
JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"

# IMPORTANT: include ALL jars under lib/ to resolve BayesianNhmmV5, etc.
CP="${JAR}:lib/*"

MODEL="data/zenodo_ref/healthy_WGS.mincg7.example.hmm_model"
IN1="out/CpgMultiMetricsStats.hg19.details.bed.gz"

OUT_DIR="out"
OUTP="${OUT_DIR}/test.prediction.bed.gz"
TMP_OUT="${OUTP}.tmp.$$"
LOG="${OUT_DIR}/step3_decode.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

# --- sanity checks ---
for f in "$JAR" "$MODEL" "$IN1"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

# avoid leaving broken gzip on failure
rm -f "$TMP_OUT"

# --- Step3: decode ---
java -Xmx10G -cp "$CP" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe \
  "$MODEL" \
  "$IN1" \
  "$TMP_OUT" \
  -decodeModeOnly \
  2> "$LOG"

# validate gzip then atomically publish
gzip -t "$TMP_OUT"
mv -f "$TMP_OUT" "$OUTP"

echo "OK: wrote $OUTP" >&2
echo "LOG: $LOG" >&2
