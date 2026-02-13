#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

# ===== Java classpath =====
JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/*"

# ===== Inputs =====
PARTS_DIR="out/parts"
# 学習に使う part 数（小さめで十分。まずは 8 くらい推奨）
TRAIN_PART_N="${1:-8}"

# ===== Outputs =====
OUT_DIR="out"
mkdir -p "$OUT_DIR"

MODEL_OUT="${OUT_DIR}/selftrained.states2.features3.hmm_model"
PRED_TRAIN="${OUT_DIR}/training.pred.gz"
LOG="${OUT_DIR}/step2_train_refit.$(date +%Y%m%d_%H%M%S).log"

# ===== FinaleMe training params =====
XMX_GB="${2:-24}"
FEATURES="${3:-3}"
STATES="${4:-2}"
ITER="${5:-20}"
MINI_DP="${6:-2}"

# ===== Checks =====
for f in "$JAR"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

# collect train parts (part000..)
TRAIN_PARTS=()
for i in $(seq -w 0 $((TRAIN_PART_N-1))); do
  f="${PARTS_DIR}/input_matrix.part${i}.tsv.gz"
  if [ ! -f "$f" ]; then
    echo "ERROR: missing train part: $f" >&2
    exit 1
  fi
  TRAIN_PARTS+=( "$f" )
done

echo "INFO: TRAIN_PART_N=$TRAIN_PART_N" >&2
echo "INFO: XMX=${XMX_GB}G FEATURES=$FEATURES STATES=$STATES ITER=$ITER MINI_DP=$MINI_DP" >&2
echo "INFO: MODEL_OUT=$MODEL_OUT" >&2
echo "INFO: PRED_TRAIN=$PRED_TRAIN" >&2
echo "INFO: LOG=$LOG" >&2

# ===== Build a single training input (concatenate parts) =====
# FinaleMe は1ファイル入力想定なので、まず学習用に結合する
TRAIN_INPUT="${OUT_DIR}/input_matrix.train.parts0-$(printf "%03d" $((TRAIN_PART_N-1))).tsv.gz"
TMP_TRAIN="${TRAIN_INPUT}.tmp.$$"
rm -f "$TMP_TRAIN"

# ヘッダは最初の part の1行目だけ採用し、以降の part はヘッダを落として結合
{
  zcat "${TRAIN_PARTS[0]}" | head -n 1
  for f in "${TRAIN_PARTS[@]}"; do
    zcat "$f" | tail -n +2
  done
} | gzip -c > "$TMP_TRAIN"

gzip -t "$TMP_TRAIN"
mv -f "$TMP_TRAIN" "$TRAIN_INPUT"

echo "OK: wrote TRAIN_INPUT=$TRAIN_INPUT" >&2
ls -lh "$TRAIN_INPUT" >&2
echo -n "INFO: TRAIN_INPUT lines=" >&2
zcat "$TRAIN_INPUT" | wc -l >&2

# ===== Train (and also write training prediction) =====
# 重要: decode-only は付けない（学習を実行する）
java -Xmx"${XMX_GB}G" -cp "$CP" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe \
  -features "$FEATURES" \
  -states "$STATES" \
  -iteration "$ITER" \
  -miniDataPoints "$MINI_DP" \
  "$MODEL_OUT" \
  "$TRAIN_INPUT" \
  "$PRED_TRAIN" \
  2> "$LOG"

echo "OK: wrote $MODEL_OUT" >&2
echo "OK: wrote $PRED_TRAIN" >&2
echo "LOG: $LOG" >&2
