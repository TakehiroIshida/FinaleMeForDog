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
#
# Examples:
#   ./second_step_from_detail.sh 50000
#   ./second_step_from_detail.sh 50000 24 2 20
#   ./second_step_from_detail.sh 50000 24 2 10   # pick more fragments (less strict)
#
# Notes:
# - FinaleMe training forces miniDataPoints >= 2 internally, so MINI_DP=1 is not meaningful.
N_FRAGS="${1:-50000}"
XMX_GB="${2:-12}"
MINI_DP="${3:-2}"
HASH_DIV="${4:-20}"   # smaller => keep more readNames, larger => keep fewer

# keep details format (header + 12 columns)
INPUT_MAT="${OUT_DIR}/input_matrix.details.hashN${N_FRAGS}.div${HASH_DIV}.tsv.gz"

# --- sanity checks ---
for f in "$JAR" "$DETAILS"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

# --- build input by "scattering" readName selection + pre-filtering to match FinaleMe ---
# Details columns:
# 1 chr, 2 start, 3 end, 4 readName, 5 FragLen, 6 Frag_strand, 7 methy_stat,
# 8 Norm_Frag_cov, 9 baseQ, 10 Offset_frag, 11 Dist_frag_end, 12 methyPrior
TMP_INPUT="${INPUT_MAT}.tmp.$$"
rm -f "$TMP_INPUT"

zcat "$DETAILS" \
| awk -v N="$N_FRAGS" -v DIV="$HASH_DIV" 'BEGIN{FS=OFS="\t"}
  NR==1{print; next}

  {
    # ---- line-level filters to align with FinaleMe.processMatrixFile() ----
    fraglen=$5+0
    baseQ=$9+0
    offset=$10+0
    prior=$12

    # same style as FinaleMe defaults (minFragLen=30, maxFragLen=500, baseQ>5, offset>=0)
    if (fraglen >= 500 || fraglen <= 30) next
    if (baseQ <= 5) next
    if (offset < 0) next

    # methyPrior NaN rows are skipped in FinaleMe, so fix them here (keep column count)
    if (prior=="NaN" || prior=="nan" || prior=="") $12=0

    rn=$4

    # ---- select up to N unique readNames, but "scatter" selection by hashing readName ----
    # readName is numeric in your data; use modulo as cheap hash.
    if (!(rn in keep)) {
      # keep roughly 1/DIV of readNames to avoid "file head bias"
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
      print "WARN: selected fragments < N_FRAGS. Consider lowering HASH_DIV (e.g., 10 or 5)." > "/dev/stderr"
    }
  }' \
| gzip -c > "$TMP_INPUT"

gzip -t "$TMP_INPUT"
mv -f "$TMP_INPUT" "$INPUT_MAT"

echo "OK: wrote $INPUT_MAT" >&2

# --- quick checks ---
echo "INFO: unique readName count (post-filter):" >&2
zcat "$INPUT_MAT" | awk 'BEGIN{FS="\t"} NR>1{u[$4]=1} END{print length(u)}' >&2

echo "INFO: top readName multiplicities (post-filter):" >&2
zcat "$INPUT_MAT" | awk 'BEGIN{FS="\t"} NR>1{c[$4]++} END{for(k in c) print c[k],k}' | sort -nr | head >&2

# --- Step2: train ---
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
