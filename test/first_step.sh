#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}/.."

JAR="target/FinaleMe-0.58-jar-with-dependencies.jar"
CP="${JAR}:lib/gatk-package-distribution-3.3.jar:lib/sis-jhdf5-batteries_included.jar:lib/java-genomics-io.jar:lib/igv.jar"

CG="data/zenodo_ref/CG_motif.hg19.common_chr.pos_only.bedgraph.gz"
DARK="data/zenodo_ref/dark.hg19.bed"
BAM="test/data/BH01.chr22.chr.bam"
REF2BIT="data/zenodo_ref/hg19.2bit"
PRIOR="data/zenodo_ref/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw"

OUT_DIR="out"
OUT1="${OUT_DIR}/CpgMultiMetricsStats.hg19.details.bed.gz"
TMP_OUT="${OUT1}.tmp.$$"
LOG="${OUT_DIR}/first_step.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUT_DIR"

for f in "$JAR" "${REF2BIT}" "${CG}" "${DARK}" "${PRIOR}" "${BAM}"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

# BAM index 明示チェック
if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
  echo "ERROR: BAM index (.bai) not found for $BAM" >&2
  exit 1
fi

rm -f "$TMP_OUT"

java -Xmx24G -cp "$CP" org.cchmc.epifluidlab.finaleme.utils.CpgMultiMetricsStats \
  "$REF2BIT" \
  "$CG" \
  "$CG" \
  "$BAM" \
  "$TMP_OUT" \
  -excludeRegions "$DARK" \
  -valueWigs methyPrior:0:"$PRIOR" \
  -wgsMode \
  2> "$LOG"

gzip -t "$TMP_OUT"
mv -f "$TMP_OUT" "$OUT1"
echo "OK: wrote $OUT1" >&2
echo "LOG: $LOG" >&2
