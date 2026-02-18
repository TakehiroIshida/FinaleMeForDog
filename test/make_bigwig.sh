#!/usr/bin/env bash
set -euo pipefail

IN="out/prediction.all.parts128.gz"   # ← 今の実ファイル
CHRSIZES="test/data/BH01.chrom.sizes"
PREFIX="BH01"

OUTDIR="out/bedgraph"
mkdir -p "$OUTDIR"

BW_IMG="$HOME/work/images/ucsc-bedgraphtobigwig.sif"

echo "INFO: using image $BW_IMG" >&2

# -------- bedGraph 作成 --------

make_bg() {
  local col="$1"
  local out="$2"

  zcat "$IN" \
    | grep -v '^#' \
    | sed 's/^chrM/MT/; s/^chr//' \
    | cut -f1-3,"$col" \
    | bedtools sort -i stdin \
    > "$out"
}

echo "INFO: generating bedGraph files..." >&2

make_bg 4 "${OUTDIR}/${PREFIX}.methy_perc_predict.bedgraph"
make_bg 5 "${OUTDIR}/${PREFIX}.methy_count_predict.bedgraph"
make_bg 6 "${OUTDIR}/${PREFIX}.total_count_predict.bedgraph"

make_bg 7 "${OUTDIR}/${PREFIX}.methy_perc_obs.bedgraph"
make_bg 8 "${OUTDIR}/${PREFIX}.methy_count_obs.bedgraph"
make_bg 9 "${OUTDIR}/${PREFIX}.total_count_obs.bedgraph"

# -------- bigWig 変換（Apptainer経由）--------

bw() {
  local in_bg="$1"
  local out_bw="$2"

  echo "INFO: converting $in_bg -> $out_bw" >&2
  apptainer exec -e "$BW_IMG" \
    bedGraphToBigWig "$in_bg" "$CHRSIZES" "$out_bw"
}

echo "INFO: converting to bigWig..." >&2

bw "${OUTDIR}/${PREFIX}.methy_perc_predict.bedgraph"  "${OUTDIR}/${PREFIX}.methy_perc_predict.bw"
bw "${OUTDIR}/${PREFIX}.methy_count_predict.bedgraph" "${OUTDIR}/${PREFIX}.methy_count_predict.bw"
bw "${OUTDIR}/${PREFIX}.total_count_predict.bedgraph" "${OUTDIR}/${PREFIX}.total_count_predict.bw"

bw "${OUTDIR}/${PREFIX}.methy_perc_obs.bedgraph"      "${OUTDIR}/${PREFIX}.methy_perc_obs.bw"
bw "${OUTDIR}/${PREFIX}.methy_count_obs.bedgraph"     "${OUTDIR}/${PREFIX}.methy_count_obs.bw"
bw "${OUTDIR}/${PREFIX}.total_count_obs.bedgraph"     "${OUTDIR}/${PREFIX}.total_count_obs.bw"

echo "OK: done." >&2
ls -lh "$OUTDIR"/*.bw >&2
