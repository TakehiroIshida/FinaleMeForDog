#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

# ===== Java classpath =====
JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/*"

# ===== Inputs =====
PARTS_DIR="out/parts"
TRAIN_PART_N="${1:-8}"     # 学習に使うpart数（まず8推奨）

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
if [ ! -f "$JAR" ]; then
  echo "ERROR: missing file: $JAR" >&2
  exit 1
fi

# collect train parts: part000..part(N-1)
TRAIN_PARTS=()
for ((i=0; i<TRAIN_PART_N; i++)); do
  part="$(printf "%03d" "$i")"
  f="${PARTS_DIR}/input_matrix.part${part}.tsv.gz"
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
# FinaleMe は1ファイル入力なので、学習用に結合ファイルを作る
LAST_PART="$(printf "%03d" $((TRAIN_PART_N-1)))"
TRAIN_INPUT="${OUT_DIR}/input_matrix.train.parts000-${LAST_PART}.tsv.gz"
TMP_TRAIN="${TRAIN_INPUT}.tmp.$$"
rm -f "$TMP_TRAIN"

# header を安全に取得（パイプで head を使わない）
HEADER="$(python - << 'PY'
import gzip
p="out/parts/input_matrix.part000.tsv.gz"
with gzip.open(p, "rt") as f:
    print(f.readline().rstrip("\n"))
PY
)"

# 結合：先頭に header、各partは2行目以降を足す（awk NR>1）
{
  printf '%s\n' "$HEADER"
  for f in "${TRAIN_PARTS[@]}"; do
    zcat "$f" | awk 'NR>1'
  done
} | gzip -c > "$TMP_TRAIN"

gzip -t "$TMP_TRAIN"
mv -f "$TMP_TRAIN" "$TRAIN_INPUT"

echo "OK: wrote TRAIN_INPUT=$TRAIN_INPUT" >&2
ls -lh "$TRAIN_INPUT" >&2
echo -n "INFO: TRAIN_INPUT lines=" >&2
zcat "$TRAIN_INPUT" | wc -l >&2

# ===== Train (and also write training prediction) =====
# stdout/stderr をどちらもログに保存
java -Xmx"${XMX_GB}G" -cp "$CP" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe \
  -features "$FEATURES" \
  -states "$STATES" \
  -iteration "$ITER" \
  -miniDataPoints "$MINI_DP" \
  "$MODEL_OUT" \
  "$TRAIN_INPUT" \
  "$PRED_TRAIN" \
  > "$LOG" 2>&1

echo "OK: wrote $MODEL_OUT" >&2
echo "OK: wrote $PRED_TRAIN" >&2
echo "LOG: $LOG" >&2
