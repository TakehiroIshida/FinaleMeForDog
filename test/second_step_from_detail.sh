#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/*"

# Step1 output（既にできている前提）
DETAILS="out/CpgMultiMetricsStats.hg19.details.bed.gz"

OUT_DIR="out"
# ここは「detailsをそのまま学習入力として渡す」ので、名前だけ input_matrix にしておく
INPUT_MAT="${OUT_DIR}/input_matrix.details.tsv.gz"
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

# --- build input matrix: keep DETAILS columns (drop header only) ---
TMP_INPUT="${INPUT_MAT}.tmp.$$"
rm -f "$TMP_INPUT"

# details はヘッダ付きなので除去。
# methyPrior の NaN が気になる場合もあるので 0 に置換（列数は維持）
zcat "$DETAILS" \
| awk 'BEGIN{FS=OFS="\t"}
  NR==1{next}
  {
    # methyPrior は末尾列にいる想定。空/NaNなら0にする（列数を変えない）
    if ($NF=="NaN" || $NF=="nan" || $NF=="") $NF=0;
    print
  }' \
| gzip -c > "$TMP_INPUT"

gzip -t "$TMP_INPUT"
mv -f "$TMP_INPUT" "$INPUT_MAT"

echo "OK: wrote $INPUT_MAT" >&2

# --- Step2: train (no -decodeModeOnly) ---
# 完走優先: miniDataPoints を 1 に落として全除外を回避
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
