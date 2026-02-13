IN="out/prediction.all.parts128.gz"
CHRSIZES="test/data/BH01.chrom.sizes"
PREFIX="BH01"

# 予測（predict）
zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,4 \
  | bedtools sort -i stdin \
  > "${PREFIX}.methy_perc_predict.bedgraph"

zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,5 \
  | bedtools sort -i stdin \
  > "${PREFIX}.methy_count_predict.bedgraph"

zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,6 \
  | bedtools sort -i stdin \
  > "${PREFIX}.total_count_predict.bedgraph"

# 観測（obs）
zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,7 \
  | bedtools sort -i stdin \
  > "${PREFIX}.methy_perc_obs.bedgraph"

zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,8 \
  | bedtools sort -i stdin \
  > "${PREFIX}.methy_count_obs.bedgraph"

zcat "$IN" \
  | grep -v '^#' \
  | sed 's/^chrM/MT/; s/^chr//' \
  | cut -f1-3,9 \
  | bedtools sort -i stdin \
  > "${PREFIX}.total_count_obs.bedgraph"

# bedGraph -> bigWig
bedGraphToBigWig "${PREFIX}.methy_perc_predict.bedgraph" "$CHRSIZES" "${PREFIX}.methy_perc_predict.bw"
bedGraphToBigWig "${PREFIX}.methy_count_predict.bedgraph" "$CHRSIZES" "${PREFIX}.methy_count_predict.bw"
bedGraphToBigWig "${PREFIX}.total_count_predict.bedgraph" "$CHRSIZES" "${PREFIX}.total_count_predict.bw"

bedGraphToBigWig "${PREFIX}.methy_perc_obs.bedgraph" "$CHRSIZES" "${PREFIX}.methy_perc_obs.bw"
bedGraphToBigWig "${PREFIX}.methy_count_obs.bedgraph" "$CHRSIZES" "${PREFIX}.methy_count_obs.bw"
bedGraphToBigWig "${PREFIX}.total_count_obs.bedgraph" "$CHRSIZES" "${PREFIX}.total_count_obs.bw"
