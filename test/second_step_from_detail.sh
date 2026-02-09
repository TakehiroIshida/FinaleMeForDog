#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/*"

DETAILS="out/CpgMultiMetricsStats.hg19.details.bed.gz"

OUT_DIR="out"
MODEL_OUT="${OUT_DIR}/selftrained.states2.features3.hmm_model"
PRED_TRAIN="${OUT_DIR}/training.pred.gz"
LOG="${OUT_DIR}/step2_train.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

# Usage:
#   ./second_step_from_detail.sh [N_FRAGS] [XMX_GB] [MINI_DP] [HASH_DIV] [MAX_PER_READ]
#
# Examples:
#   ./second_step_from_detail.sh 5000 12 2 5 20
#   ./second_step_from_detail.sh 10000 24 2 5 30
N_FRAGS="${1:-5000}"
XMX_GB="${2:-12}"
MINI_DP="${3:-2}"
HASH_DIV="${4:-5}"
MAX_PER_READ="${5:-20}"   # 1フラグメントから最大何行出すか（OOM対策）

INPUT_MAT="${OUT_DIR}/input_matrix.details.hashN${N_FRAGS}.div${HASH_DIV}.max${MAX_PER_READ}.tsv.gz"

for f in "$JAR" "$DETAILS"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

TMP_INPUT="${INPUT_MAT}.tmp.$$"
rm -f "$TMP_INPUT"

zcat "$DETAILS" \
| awk -v N="$N_FRAGS" -v DIV="$HASH_DIV" -v MAXR="$MAX_PER_READ" 'BEGIN{FS=OFS="\t"}
  NR==1{print; next}

  # 列ズレ行を排除
  NF!=12{next}

  {
    fraglen=$5+0
    baseQ=$9+0
    offset=$10+0
    prior=$12

    # FinaleMe 側の行フィルタに寄せる :contentReference[oaicite:1]{index=1}
    if (fraglen >= 500 || fraglen <= 30) next
    if (baseQ <= 5) next
    if (offset < 0) next
    if (prior=="NaN" || prior=="nan" || prior=="") $12=0

    rn=$4

    # readName を散らして採用（先頭偏り回避）
    if (!(rn in keep)) {
      if ((rn % DIV) != 0) next
      if (k >= N) next
      keep[rn]=1
      k++
    }

    # OOM対策：1フラグメントあたりの出力行数を制限
    if (++seen[rn] > MAXR) next

    print
  }
  END{
    print "INFO: selected fragments =", k > "/dev/stderr"
    total=0
    for (x in seen) total += seen[x]
    print "INFO: total lines written =", total > "/dev/stderr"
  }' \
| gzip -c > "$TMP_INPUT"

gzip -t "$TMP_INPUT"
mv -f "$TMP_INPUT" "$INPUT_MAT"
echo "OK: wrote $INPUT_MAT" >&2

java -Xmx"${XMX_GB}G" -cp "$CP" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe \
  -features 3 \
  -states 2 \
  -iteration 20 \
  -miniDataPoints "$MINI_DP" \
  "$MODEL_OUT" \
  "$INPUT_MAT" \
  "$PRED_TRAIN" \
  2> "$LOG"

echo "OK: wrote $MODEL_OUT" >&2
echo "OK: wrote $PRED_TRAIN" >&2
echo "LOG: $LOG" >&2
