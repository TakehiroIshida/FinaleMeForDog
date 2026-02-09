#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/*"

# Step1 output（既にできている前提）
DETAILS="out/CpgMultiMetricsStats.hg19.details.bed.gz"

OUT_DIR="out"
INPUT_MAT="${OUT_DIR}/input_matrix.features3.tsv.gz"
MODEL_OUT="${OUT_DIR}/selftrained.states2.features3.hmm_model"
PRED_TRAIN="${OUT_DIR}/training.pred.gz"
LOG="${OUT_DIR}/step2_train.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

# --- sanity checks ---
for f in "$JAR" "$DETAILS"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

# --- build input matrix (chr, start, end, fraglen, cov, prior) ---
TMP_INPUT="${INPUT_MAT}.tmp.$$"
rm -f "$TMP_INPUT"

zcat "$DETAILS" \
| awk 'BEGIN{OFS="\t"}
  NR==1{next}
  {
    chr=$1; start=$2; end=$3;
    fraglen=$5+0;
    cov=$8+0;
    prior=$12;
    if (prior=="NaN" || prior=="nan" || prior=="") prior=0;
    print chr, start, end, fraglen, cov, prior
  }' \
| gzip -c > "$TMP_INPUT"

gzip -t "$TMP_INPUT"
mv -f "$TMP_INPUT" "$INPUT_MAT"

echo "OK: wrote $INPUT_MAT" >&2

# --- Step2: train (no -decodeModeOnly) ---
# まず完走優先で iteration は20、statesは2、featuresは3
java -Xmx12G -cp "$CP" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe \
  -features 3 \
  -states 2 \
  -iteration 20 \
  "$MODEL_OUT" \
  "$INPUT_MAT" \
  "$PRED_TRAIN" \
  2> "$LOG"

echo "OK: wrote $MODEL_OUT" >&2
echo "OK: wrote $PRED_TRAIN" >&2
echo "LOG: $LOG" >&2
