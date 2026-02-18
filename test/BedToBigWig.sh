#!/usr/bin/env bash
set -euo pipefail

IN="out/prediction.all.parts.gz"
CHRSIZES="test/data/BH01.chrom.sizes"
PREFIX="BH01"

# 出力先（FinaleMeForDog/out 配下）
OUTDIR="out/bedgraph"
mkdir -p "$OUTDIR"

# 予測（predict）
zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,4 \
  | bedtools sort -i stdin \
  > "${OUTDIR}/${PREFIX}.methy_perc_predict.bedgraph"

zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,5 \
  | bedtools sort -i stdin \
  > "${OUTDIR}/${PREFIX}.methy_count_predict.bedgraph"

zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,6 \
  | bedtools sort -i stdin \
  > "${OUTDIR}/${PREFIX}.total_count_predict.bedgraph"

# 観測（obs）
zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,7 \
  | bedtools sort -i stdin \
  > "${OUTDIR}/${PREFIX}.methy_perc_obs.bedgraph"

zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,8 \
  | bedtools sort -i stdin \
  > "${OUTDIR}/${PREFIX}.methy_count_obs.bedgraph"

zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,9 \
  | bedtools sort -i stdin \
  > "${OUTDIR}/${PREFIX}.total_count_obs.bedgraph"

# bedGraph -> bigWig（同じ out/bedgraph に出す）
bedGraphToBigWig "${OUTDIR}/${PREFIX}.methy_perc_predict.bedgraph" "$CHRSIZES" "${OUTDIR}/${PREFIX}.methy_perc_predict.bw"
bedGraphToBigWig "${OUTDIR}/${PREFIX}.methy_count_predict.bedgraph" "$CHRSIZES" "${OUTDIR}/${PREFIX}.methy_count_predict.bw"
bedGraphToBigWig "${OUTDIR}/${PREFIX}.total_count_predict.bedgraph" "$CHRSIZES" "${OUTDIR}/${PREFIX}.total_count_predict.bw"

bedGraphToBigWig "${OUTDIR}/${PREFIX}.methy_perc_obs.bedgraph" "$CHRSIZES" "${OUTDIR}/${PREFIX}.methy_perc_obs.bw"
bedGraphToBigWig "${OUTDIR}/${PREFIX}.methy_count_obs.bedgraph" "$CHRSIZES" "${OUTDIR}/${PREFIX}.methy_count_obs.bw"
bedGraphToBigWig "${OUTDIR}/${PREFIX}.total_count_obs.bedgraph" "$CHRSIZES" "${OUTDIR}/${PREFIX}.total_count_obs.bw"
