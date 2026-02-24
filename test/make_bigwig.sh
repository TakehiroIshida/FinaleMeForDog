#!/usr/bin/env bash
set -euo pipefail

IN="out/prediction.all.G1000.low_cov.gz"
CHRSIZES="test/data/BH01.chrom.sizes"
PREFIX="BH01"

OUTDIR="out/bedgraph"
mkdir -p "$OUTDIR"

BW_IMG="$HOME/work/images/ucsc-bedgraphtobigwig.sif"

[ -f "$IN" ] || { echo "ERROR: missing IN=$IN" >&2; exit 1; }
[ -f "$CHRSIZES" ] || { echo "ERROR: missing CHRSIZES=$CHRSIZES" >&2; exit 1; }
[ -f "$BW_IMG" ] || { echo "ERROR: missing BW_IMG=$BW_IMG" >&2; exit 1; }

# BigWig (Apptainer)
bw() {
  local in_bg="$1"
  local out_bw="$2"
  apptainer exec -e "$BW_IMG" bedGraphToBigWig "$in_bg" "$CHRSIZES" "$out_bw"
}

# ---- 1) カウント系を座標ごとに合計して bedGraph を作る ----
# 出力:
#  - methy_count_predict.bedgraph (col5 sum)
#  - total_count_predict.bedgraph (col6 sum)
#  - methy_count_obs.bedgraph     (col8 sum)
#  - total_count_obs.bedgraph     (col9 sum)
make_sum_bg() {
  local col="$1"
  local out="$2"

  zcat "$IN" \
    | grep -v '^#' \
    | sed 's/^chrM/MT/; s/^chr//' \
    | cut -f1-3,"$col" \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    | awk 'BEGIN{OFS="\t"}
        NR==1{c=$1;s=$2;e=$3;sum=$4;next}
        ($1==c && $2==s && $3==e){sum+=$4;next}
        {print c,s,e,sum; c=$1;s=$2;e=$3;sum=$4}
        END{if(NR>0) print c,s,e,sum}
      ' \
    > "$out"
}

echo "INFO: summing counts..." >&2
make_sum_bg 5 "${OUTDIR}/${PREFIX}.methy_count_predict.bedgraph"
make_sum_bg 6 "${OUTDIR}/${PREFIX}.total_count_predict.bedgraph"
make_sum_bg 8 "${OUTDIR}/${PREFIX}.methy_count_obs.bedgraph"
make_sum_bg 9 "${OUTDIR}/${PREFIX}.total_count_obs.bedgraph"

# ---- 2) 割合 bedGraph を count から再計算（加重平均相当）----
# methy_perc = 100 * methy_count / total_count（total_count=0 は 0）
make_perc_from_counts() {
  local methy_bg="$1"
  local total_bg="$2"
  local out_perc="$3"

  join -t $'\t' -1 1 -2 1 \
    <(awk 'BEGIN{OFS="\t"}{print $1":"$2":"$3, $4}' "$methy_bg") \
    <(awk 'BEGIN{OFS="\t"}{print $1":"$2":"$3, $4}' "$total_bg") \
  | awk -F'\t' 'BEGIN{OFS="\t"}
      {
        split($1, a, ":");
        chr=a[1]; start=a[2]; end=a[3];
        m=$2; t=$3;
        perc = (t>0 ? (100.0*m/t) : 0.0);
        print chr, start, end, perc
      }
    ' \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
  > "$out_perc"
}

echo "INFO: computing perc from counts..." >&2
make_perc_from_counts \
  "${OUTDIR}/${PREFIX}.methy_count_predict.bedgraph" \
  "${OUTDIR}/${PREFIX}.total_count_predict.bedgraph" \
  "${OUTDIR}/${PREFIX}.methy_perc_predict.bedgraph"

make_perc_from_counts \
  "${OUTDIR}/${PREFIX}.methy_count_obs.bedgraph" \
  "${OUTDIR}/${PREFIX}.total_count_obs.bedgraph" \
  "${OUTDIR}/${PREFIX}.methy_perc_obs.bedgraph"

# ---- 3) BigWig 変換 ----
echo "INFO: converting to bigWig..." >&2

bw "${OUTDIR}/${PREFIX}.methy_perc_predict.bedgraph"  "${OUTDIR}/${PREFIX}.methy_perc_predict.bw"
bw "${OUTDIR}/${PREFIX}.methy_count_predict.bedgraph" "${OUTDIR}/${PREFIX}.methy_count_predict.bw"
bw "${OUTDIR}/${PREFIX}.total_count_predict.bedgraph" "${OUTDIR}/${PREFIX}.total_count_predict.bw"

bw "${OUTDIR}/${PREFIX}.methy_perc_obs.bedgraph"      "${OUTDIR}/${PREFIX}.methy_perc_obs.bw"
bw "${OUTDIR}/${PREFIX}.methy_count_obs.bedgraph"     "${OUTDIR}/${PREFIX}.methy_count_obs.bw"
bw "${OUTDIR}/${PREFIX}.total_count_obs.bedgraph"     "${OUTDIR}/${PREFIX}.total_count_obs.bw"

echo "OK: wrote outputs in $OUTDIR" >&2
ls -lh "$OUTDIR"/${PREFIX}.*.{bedgraph,bw} 2>/dev/null | sed 's/^/  /' >&2 || true
