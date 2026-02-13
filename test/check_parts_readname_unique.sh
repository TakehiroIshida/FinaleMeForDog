#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

OUT_DIR="out/parts"
PATTERN="${OUT_DIR}/input_matrix.part*.tsv.gz"

TMP_DIR="out/tmp_check_parts"
mkdir -p "$TMP_DIR"

# どのくらいメモリ/時間を使うかを抑えるため、まずは “各partのreadName集合のハッシュ” で粗チェックし、
# その後に厳密チェック（sort/uniq）を行う。

echo "INFO: checking part files: $PATTERN" >&2

PART_FILES=( $PATTERN )
if [ "${#PART_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no part files found: $PATTERN" >&2
  exit 1
fi

echo "INFO: found ${#PART_FILES[@]} part files" >&2

# --- 1) 各partのreadNameを抽出してソート（ヘッダ除外） ---
# 4列目がreadName
# 出力: out/tmp_check_parts/partXXX.reads.sorted.txt
for f in "${PART_FILES[@]}"; do
  base="$(basename "$f")"                     # input_matrix.part000.tsv.gz
  part="${base#input_matrix.part}"            # 000.tsv.gz
  part="${part%%.*}"                          # 000

  out_reads="${TMP_DIR}/part${part}.reads.sorted.txt"

  echo "INFO: extracting readNames part=${part}" >&2
  zcat "$f" \
    | awk 'BEGIN{FS="\t"} NR==1{next} {print $4}' \
    | sort -S 1G \
    | uniq \
    > "$out_reads"
done

# --- 2) 全partのreadNameを結合し、重複がないか確認 ---
ALL_READS="${TMP_DIR}/all.reads.sorted.txt"
DUP_READS="${TMP_DIR}/dup.reads.txt"

echo "INFO: merging and checking duplicates..." >&2

cat "${TMP_DIR}"/part*.reads.sorted.txt \
  | sort -S 2G \
  > "$ALL_READS"

# 重複readName（複数partに存在する）を抽出
uniq -d "$ALL_READS" > "$DUP_READS"

DUP_N=$(wc -l < "$DUP_READS" | tr -d ' ')
if [ "$DUP_N" -ne 0 ]; then
  echo "ERROR: duplicated readName across parts detected: $DUP_N" >&2
  echo "INFO: examples (first 20):" >&2
  head -n 20 "$DUP_READS" >&2
  echo "INFO: full list: $DUP_READS" >&2
  exit 1
fi

echo "OK: no duplicated readName across parts" >&2

# --- 3) 後始末（必要なら消す） ---
echo "INFO: tmp outputs under: $TMP_DIR" >&2
echo "INFO: if you want to remove: rm -r $TMP_DIR" >&2
