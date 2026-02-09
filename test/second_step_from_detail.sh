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
#   ./second_step_from_detail.sh [N_FRAGS] [XMX_GB] [MINI_DP] [HASH_DIV]
N_FRAGS="${1:-50000}"
XMX_GB="${2:-12}"
MINI_DP="${3:-2}"        # trainingでは2未満は2に引き上げられる :contentReference[oaicite:2]{index=2}
HASH_DIV="${4:-10}"      # 小さくするとreadNameを多めに拾う

INPUT_MAT="${OUT_DIR}/input_matrix.details.hashN${N_FRAGS}.div${HASH_DIV}.tsv.gz"

for f in "$JAR" "$DETAILS"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

TMP_INPUT="${INPUT_MAT}.tmp.$$"
rm -f "$TMP_INPUT"

# details header:
# chr start end readName FragLen Frag_strand methy_stat Norm_Frag_cov baseQ Offset_frag Dist_frag_end methyPrior
zcat "$DETAILS" \
| awk -v N="$N_FRAGS" -v DIV="$HASH_DIV" 'BEGIN{FS=OFS="\t"}
  NR==1{print; next}

  # 列ズレ行を排除（重要）
  NF!=12{next}

  {
    fraglen=$5+0
    baseQ=$9+0
    offset=$10+0
    prior=$12

    # FinaleMe.processMatrixFile() の行フィルタに合わせる :contentReference[oaicite:3]{index=3}
    if (fraglen >= 500 || fraglen <= 30) next
    if (baseQ <= 5) next
    if (offset < 0) next

    # FinaleMeはpriorがNaNだとその行を捨てるので0に補正 :contentReference[oaicite:4]{index=4}
    if (prior=="NaN" || prior=="nan" || prior=="") $12=0

    rn=$4

    # readNameを先頭偏りなく拾う（readNameが数値の前提）
    if (!(rn in keep)) {
      if ((rn % DIV) != 0) next
      if (k >= N) next
      keep[rn]=1
      k++
    }

    print
  }
  END{
    print "INFO: selected fragments =", k > "/dev/stderr"
    if (k < N) {
      print "WARN: selected fragments < N_FRAGS. Lower HASH_DIV (e.g., 5) to pick more." > "/dev/stderr"
    }
  }' \
| gzip -c > "$TMP_INPUT"

gzip -t "$TMP_INPUT"
mv -f "$TMP_INPUT" "$INPUT_MAT"
echo "OK: wrote $INPUT_MAT" >&2

echo "INFO: unique readName count (post-filter):" >&2
zcat "$INPUT_MAT" | awk 'BEGIN{FS="\t"} NR>1{u[$4]=1} END{print length(u)}' >&2

echo "INFO: top readName multiplicities (post-filter):" >&2
zcat "$INPUT_MAT" | awk 'BEGIN{FS="\t"} NR>1{c[$4]++} END{for(k in c) print c[k],k}' | sort -nr | head >&2

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
