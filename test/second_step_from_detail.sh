#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/*"

# Step1 output（既にできている前提）
DETAILS="out/CpgMultiMetricsStats.hg19.details.bed.gz"

OUT_DIR="out"
MODEL_OUT="${OUT_DIR}/selftrained.states2.features3.hmm_model"
PRED_TRAIN="${OUT_DIR}/training.pred.gz"
LOG="${OUT_DIR}/step2_train.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

# -------- parameters --------
# Usage:
#   ./second_step_from_detail.sh [N_FRAGS] [XMX_GB]
# Examples:
#   ./second_step_from_detail.sh 50000
#   ./second_step_from_detail.sh 200000 24
N_FRAGS="${1:-50000}"
XMX_GB="${2:-12}"

# readName単位で抽出した学習用入力（details形式を維持）
INPUT_MAT="${OUT_DIR}/input_matrix.details.first${N_FRAGS}frags.tsv.gz"

# --- sanity checks ---
for f in "$JAR" "$DETAILS"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

# --- build input matrix by sampling first N_FRAGS fragments (keep series by readName) ---
TMP_INPUT="${INPUT_MAT}.tmp.$$"
rm -f "$TMP_INPUT"

zcat "$DETAILS" \
| awk -v N="$N_FRAGS" 'BEGIN{FS=OFS="\t"}
  NR==1{next}
  {
    rn=$4
    if (!(rn in keep)) {
      if (k >= N) next
      keep[rn]=1
      k++
    }
    # methyPrior is last column; keep column count, replace missing/NaN with 0
    if ($NF=="NaN" || $NF=="nan" || $NF=="") $NF=0
    print
  }
  END{
    # progress info to stderr
    print "INFO: selected fragments =", k > "/dev/stderr"
  }' \
| gzip -c > "$TMP_INPUT"

gzip -t "$TMP_INPUT"
mv -f "$TMP_INPUT" "$INPUT_MAT"

echo "OK: wrote $INPUT_MAT" >&2

# --- quick checks (optional but cheap) ---
echo "INFO: unique readName count (<= N_FRAGS):" >&2
zcat "$INPUT_MAT" | cut -f4 | sort -u | wc -l >&2

echo "INFO: top readName multiplicities (should show >=2 if series exist):" >&2
zcat "$INPUT_MAT" | cut -f4 | sort | uniq -c | sort -nr | head >&2

# --- Step2: train (no -decodeModeOnly) ---
# 完走優先: miniDataPoints を 1 に落として全除外を回避
java -Xmx"${XMX_GB}G" -cp "$CP" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe \
  -features 3 \
  -states 2 \
  -iteration 20 \
  -miniDataPoints 1 \
  "$MODEL_OUT" \
  "$INPUT_MAT" \
  "$PRED_TRAIN" \
  2> "$LOG"

echo "OK: wrote $MODEL_OUT" >&2
echo "OK: wrote $PRED_TRAIN" >&2
echo "LOG: $LOG" >&2
