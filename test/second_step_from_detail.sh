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

# --- build input matrix (chr, start, end, readName, fraglen, cov, distToCenter) ---
# details header:
# chr start end readName FragLen Frag_strand methy_stat Norm_Frag_cov baseQ Offset_frag Dist_frag_end methyPrior
TMP_INPUT="${INPUT_MAT}.tmp.$$"
rm -f "$TMP_INPUT"

zcat "$DETAILS" \
| awk 'BEGIN{FS=OFS="\t"}
  NR==1{next}
  {
    chr=$1; start=$2; end=$3; read=$4;
    fraglen=$5+0;
    cov=$8+0;
    offset=$10+0;

    # minimal guards
    if (read=="" || fraglen<=0) next;

    # distance to fragment center (integer bp):
    # center = floor(fraglen/2)
    center = int(fraglen/2);
    dist = offset - center;
    if (dist < 0) dist = -dist;

    print chr, start, end, read, fraglen, cov, dist
  }' \
| gzip -c > "$TMP_INPUT"

gzip -t "$TMP_INPUT"
mv -f "$TMP_INPUT" "$INPUT_MAT"

echo "OK: wrote $INPUT_MAT" >&2

# --- Step2: train (no -decodeModeOnly) ---
# 完走優先: iteration=20, states=2, features=3
# 既定のフラグメント内最小CpG数フィルタで全除外になるのを避けるため miniDataPoints を下げる
java -Xmx12G -cp "$CP" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe \
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
