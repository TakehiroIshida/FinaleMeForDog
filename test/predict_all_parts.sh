#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

# ===== Java classpath =====
JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/*"

# ===== Inputs =====
MODEL="out/selftrained.states2.features3.hmm_model"
PARTS_DIR="out/parts"
PART_GLOB="${PARTS_DIR}/input_matrix.part*.tsv.gz"

# ===== Outputs =====
PRED_DIR="out/pred_parts"
LOG_DIR="out/logs_pred_parts"
mkdir -p "$PRED_DIR" "$LOG_DIR"

# ===== Runtime =====
# まずは 12GB くらいから（partは小さいので通常十分）
XMX_GB="${1:-12}"

# ===== Checks =====
for f in "$JAR" "$MODEL"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

PART_FILES=( $PART_GLOB )
if [ "${#PART_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no part files found: $PART_GLOB" >&2
  exit 1
fi

echo "INFO: MODEL=$MODEL" >&2
echo "INFO: XMX=${XMX_GB}G" >&2
echo "INFO: parts=${#PART_FILES[@]}" >&2
echo "INFO: PRED_DIR=$PRED_DIR" >&2
echo "INFO: LOG_DIR=$LOG_DIR" >&2

# ===== Predict each part =====
fail=0
done_n=0
total_n="${#PART_FILES[@]}"

for in_part in "${PART_FILES[@]}"; do
  base="$(basename "$in_part")"               # input_matrix.part000.tsv.gz
  part="${base#input_matrix.part}"            # 000.tsv.gz
  part="${part%%.*}"                          # 000

  out_pred="${PRED_DIR}/prediction.part${part}.gz"
  log="${LOG_DIR}/predict.part${part}.$(date +%Y%m%d_%H%M%S).log"

  # 再実行しやすいように、既に出力があればスキップ
  if [ -f "$out_pred" ]; then
    echo "SKIP: part=$part (exists) $out_pred" >&2
    done_n=$((done_n+1))
    continue
  fi

  echo "RUN: [$((done_n+1))/$total_n] part=$part" >&2

  # decode-only（学習なし）
  if ! java -Xmx"${XMX_GB}G" -cp "$CP" org.cchmc.epifluidlab.finaleme.hmm.FinaleMe \
      -decodeModeOnly \
      "$MODEL" \
      "$in_part" \
      "$out_pred" \
      > "$log" 2>&1; then
    echo "ERROR: part=$part failed. log=$log" >&2
    fail=1
    break
  fi

  # gzip整合性
  if ! gzip -t "$out_pred" 2>/dev/null; then
    echo "ERROR: gzip test failed: $out_pred (log=$log)" >&2
    fail=1
    break
  fi

  done_n=$((done_n+1))
done

if [ "$fail" -ne 0 ]; then
  echo "FAILED" >&2
  exit 1
fi

echo "OK: all parts predicted" >&2

# ===== Quick summary =====
pred_count=$(ls -1 "$PRED_DIR"/prediction.part*.gz 2>/dev/null | wc -l | tr -d " ")
echo "INFO: predicted files count=$pred_count" >&2

first_pred="$(ls -1 "$PRED_DIR"/prediction.part*.gz | head -n 1)"
echo "INFO: sample: $first_pred" >&2
zcat "$first_pred" | head -n 5 >&2
