#!/bin/bash
# ===================================================================
# ==         PIPELINE PHÂN TÍCH DỮ LIỆU BRCA (DeepVariant)        ==
# ===================================================================

set -euo pipefail

# ========== CẤU HÌNH CHUNG ==========
SAMPLE_NAME="${1:?Vui lòng truyền SAMPLE_NAME}"
THREADS="${THREADS:-8}"
CLEANUP="${CLEANUP:-true}"
TARGET_PADDING="${TARGET_PADDING:-0}"  # đệm vùng BED (bp), 0 = tắt
TIMESTAMP() { date '+%Y-%m-%d %H:%M:%S'; }

# Conda envs
ENV_BRCA="BRCA"
ENV_GATK="GATK"
ENV_MQC="MQC"

# Dự án
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
LOG="${SAMPLE_DIR}/${SAMPLE_NAME}_pipeline.log"
mkdir -p "${SAMPLE_DIR}"
# exec > >(tee -i "$LOG") 2>&1

# ========== THAM CHIẾU & VÙNG MỤC TIÊU ==========
REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Known sites cho GATK BQSR
KNOWN_SNP="${PROJECT_DIR}/reference/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDEL="${PROJECT_DIR}/reference/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS_1000G_INDEL="${PROJECT_DIR}/reference/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# ========== TOOL PATHS ==========
FASTQC_BIN="fastqc"
TRIMMOMATIC_BIN="trimmomatic"
BWA_BIN="bwa"
SAMTOOLS_BIN="samtools"
GATK_BIN="gatk"
MOSDEPTH_BIN="mosdepth"
MULTIQC_BIN="multiqc"
BCFTOOLS_BIN="bcftools"
BEDTOOLS_BIN="bedtools"

# JAR / cấu hình
PICARD_JAR="/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar"
SNPEFF_HOME="${PROJECT_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
SNPEFF_DB="GRCh38.86"
SNPEFF_DATA_DIR="${SNPEFF_HOME}/data"

# Java opts
JAVA_OPTS_SNPEFF="-Xmx6g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_SNPSIFT="-Xmx4g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_PICARD="-Xmx4g -Djava.awt.headless=true -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"

# ========== EXTERNAL RESOURCES ==========
GNOMAD_VCF="${PROJECT_DIR}/reference/resources/gnomad.v4.1.panel.merged.vcf.gz"
CLINVAR_VCF="${PROJECT_DIR}/reference/resources/clinvar_20250810.vcf.gz"
THOUSANDG_VCF="${PROJECT_DIR}/reference/resources/1000g.panel.merged.vcf.gz"

# ========== DeepVariant ==========
RUN_DEEPVARIANT_BIN="${RUN_DEEPVARIANT_BIN:-run_deepvariant}"   # nếu cài native/conda
DV_DOCKER_IMAGE="${DV_DOCKER_IMAGE:-google/deepvariant:1.9.0}"  # Docker image
DV_USE_DOCKER="${DV_USE_DOCKER:-auto}"  # auto|true|false

# Đầu ra DeepVariant (DÙNG THƯ MỤC RIÊNG)
DEEP_DIR="${SAMPLE_DIR}/deepvariant"
DV_VCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.vcf.gz"
DV_GVCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.g.vcf.gz"
PADDED_BED="${DEEP_DIR}/${SAMPLE_NAME}.targets.pad${TARGET_PADDING}.bed"

# ========== DỮ LIỆU ĐẦU VÀO ==========
RAW_DIR="${PROJECT_DIR}/raw_data"
READ1="${RAW_DIR}/${SAMPLE_NAME}_1.fastq.gz"
READ2="${RAW_DIR}/${SAMPLE_NAME}_2.fastq.gz"
ADAPTER_FILE="/home/shmily/miniconda/envs/BRCA/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

# ========== THƯ MỤC KẾT QUẢ ==========
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
FASTQC_RAW_DIR="${SAMPLE_DIR}/fastqc_raw"
FASTQC_TRIM_DIR="${SAMPLE_DIR}/fastqc_trimmed"
BWA_DIR="${SAMPLE_DIR}/Bwa_alignments"
RECAL_DIR="${SAMPLE_DIR}/recal"
ANN_DIR="${SAMPLE_DIR}/snpeff"
COVERAGE_DIR="${SAMPLE_DIR}/coverage"
MULTIQC_DIR="${SAMPLE_DIR}/multiqc_report"

mkdir -p "${TRIM_DIR}" "${FASTQC_RAW_DIR}" "${FASTQC_TRIM_DIR}" "${BWA_DIR}" \
         "${RECAL_DIR}" "${DEEP_DIR}" "${ANN_DIR}" "${COVERAGE_DIR}" "${MULTIQC_DIR}"

# ========== FILE TRUNG GIAN ==========
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"
UNPAIRED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_unpaired.fastq.gz"
UNPAIRED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_unpaired.fastq.gz"

BAM="${BWA_DIR}/${SAMPLE_NAME}_aligned.bam"
SORTED_BAM="${BWA_DIR}/${SAMPLE_NAME}_sorted.bam"
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"
BAM_DEDUP_BAI="${BWA_DIR}/${SAMPLE_NAME}_dedup.bai"
DEDUP_METRICS="${BWA_DIR}/${SAMPLE_NAME}_dedup_metrics.txt"
STATS_FILE="${BWA_DIR}/${SAMPLE_NAME}_samtools_stats.txt"
FLAGSTAT_FILE="${BWA_DIR}/${SAMPLE_NAME}_samtools_flagstat.txt"

RECAL_DATA_TABLE="${RECAL_DIR}/${SAMPLE_NAME}_recal_data.table"
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

# ========== OUTPUT ANNOTATION (gắn nhãn _dv) ==========
SNPEFF_STATS="${ANN_DIR}/${SAMPLE_NAME}_dv"
ANN_VCF_SNPEFF="${ANN_DIR}/${SAMPLE_NAME}_dv_snpeff.vcf.gz"
ANN_VCF_GNOMAD="${ANN_DIR}/${SAMPLE_NAME}_dv_gnomad.vcf.gz"
ANN_VCF_CLINVAR="${ANN_DIR}/${SAMPLE_NAME}_dv_clinvar.vcf.gz"
ANN_VCF_FINAL="${ANN_DIR}/${SAMPLE_NAME}_dv_final_annotated.vcf.gz"

COVERAGE_PREFIX="${COVERAGE_DIR}/${SAMPLE_NAME}"
MULTIQC_FILENAME="${SAMPLE_NAME}_report.html"

# --- KÍCH HOẠT CONDA ---
eval "$(conda shell.bash hook)"


# GATK sequence dictionary
DICT="${REF%.*}.dict"
[ -s "${DICT}" ] || {
  echo "ℹ️ Tạo ${DICT}"
  ${GATK_BIN} CreateSequenceDictionary -R "${REF}" -O "${DICT}" || \
  java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" CreateSequenceDictionary R="${REF}" O="${DICT}"
}

# ===================================================================
# 1) QC & TRIMMING
# ===================================================================
echo "---=== BƯỚC 1: QC & Trimming ===---"
conda activate "${ENV_BRCA}"

${FASTQC_BIN} --threads ${THREADS} -o "${FASTQC_RAW_DIR}" "${READ1}" "${READ2}"

${TRIMMOMATIC_BIN} PE -threads ${THREADS} -phred33 \
  "${READ1}" "${READ2}" \
  "${TRIMMED_R1}" "${UNPAIRED_R1}" \
  "${TRIMMED_R2}" "${UNPAIRED_R2}" \
  ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

${FASTQC_BIN} --threads ${THREADS} -o "${FASTQC_TRIM_DIR}" "${TRIMMED_R1}" "${TRIMMED_R2}"

# ===================================================================
# 2) ALIGN + DEDUP
# ===================================================================
echo "---=== BƯỚC 2: Gióng hàng, Sắp xếp, Đánh dấu trùng lặp ===---"

${BWA_BIN} mem -Y -K 100000000 -t ${THREADS} \
  -R "@RG\tID:${SAMPLE_NAME}_RG\tSM:${SAMPLE_NAME}\tPL:Illumina\tLB:Lib1\tPU:${SAMPLE_NAME}_Illumina_Lib1" \
  "${REF}" "${TRIMMED_R1}" "${TRIMMED_R2}" \
  | ${SAMTOOLS_BIN} view -b -o "${BAM}" -

java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" SortSam \
  I="${BAM}" O="${SORTED_BAM}" SORT_ORDER=coordinate \
  VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000

java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" MarkDuplicates \
  I="${SORTED_BAM}" O="${BAM_DEDUP}" M="${DEDUP_METRICS}" CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000

${SAMTOOLS_BIN} stats "${BAM_DEDUP}" > "${STATS_FILE}"
${SAMTOOLS_BIN} flagstat "${BAM_DEDUP}" > "${FLAGSTAT_FILE}"

# ===================================================================
# 3) BQSR
# ===================================================================
echo "---=== BƯỚC 3: BQSR ===---"
conda activate "${ENV_GATK}"

${GATK_BIN} BaseRecalibrator \
  -I "${BAM_DEDUP}" -R "${REF}" \
  --known-sites "${KNOWN_SNP}" \
  --known-sites "${KNOWN_INDEL}" \
  --known-sites "${MILLS_1000G_INDEL}" \
  -O "${RECAL_DATA_TABLE}"

${GATK_BIN} ApplyBQSR \
  -R "${REF}" -I "${BAM_DEDUP}" \
  --bqsr-recal-file "${RECAL_DATA_TABLE}" \
  -O "${RECAL_BAM}"

${SAMTOOLS_BIN} index "${RECAL_BAM}"

# ===================================================================
# 4) DEEPVARIANT (TARGETED)
# ===================================================================
echo "---=== BƯỚC 4: DeepVariant (WES, regions=BED) ===---"
mkdir -p "${DEEP_DIR}"

REGIONS_BED="${TARGET_BED}"
if [ "${TARGET_PADDING}" -gt 0 ]; then
  if command -v "${BEDTOOLS_BIN}" >/dev/null 2>&1; then
    cut -f1,2 "${REF}.fai" > "${DEEP_DIR}/genome.sizes"
    ${BEDTOOLS_BIN} slop -b "${TARGET_PADDING}" -i "${TARGET_BED}" -g "${DEEP_DIR}/genome.sizes" \
      | ${BEDTOOLS_BIN} merge -i - > "${PADDED_BED}"
    REGIONS_BED="${PADDED_BED}"
  else
    echo "⚠️ Không có bedtools, dùng BED gốc."
  fi
fi

# Quyết định dùng native hay Docker
use_docker=false
if command -v "${RUN_DEEPVARIANT_BIN}" >/dev/null 2>&1; then
  use_docker=false
elif [ "${DV_USE_DOCKER}" = "true" ] || { [ "${DV_USE_DOCKER}" = "auto" ] && command -v docker >/dev/null 2>&1; }; then
  use_docker=true
else
  echo "❌ Không tìm thấy 'run_deepvariant' và Docker chưa sẵn sàng."
  exit 1
fi

echo "[4.1] Chạy DeepVariant..."
if [ "${use_docker}" = false ]; then
  "${RUN_DEEPVARIANT_BIN}" \
    --model_type=WES \
    --ref="${REF}" \
    --reads="${RECAL_BAM}" \
    --regions="${REGIONS_BED}" \
    --output_vcf="${DV_VCF}" \
    --output_gvcf="${DV_GVCF}" \
    --num_shards="${THREADS}" \
    --vcf_stats_report=true
else
  docker run --rm \
    -u "$(id -u)":"$(id -g)" \
    -v "${PROJECT_DIR}:${PROJECT_DIR}" \
    "${DV_DOCKER_IMAGE}" \
    /opt/deepvariant/bin/run_deepvariant \
      --model_type=WES \
      --ref="${REF}" \
      --reads="${RECAL_BAM}" \
      --regions="${REGIONS_BED}" \
      --output_vcf="${DV_VCF}" \
      --output_gvcf="${DV_GVCF}" \
      --num_shards="${THREADS}" \
      --vcf_stats_report=true
fi

tabix -f -p vcf "${DV_VCF}" || true
tabix -f -p vcf "${DV_GVCF}" || true
echo "✅ DeepVariant xong: ${DV_VCF}"

# ===================================================================
# 4.5) ANNOTATE (SnpEff + SnpSift)
# ===================================================================
echo "---=== BƯỚC 4.5: Annotate ===---"

# Đảm bảo DB SnpEff
if [ ! -d "${SNPEFF_DATA_DIR}/${SNPEFF_DB}" ]; then
  echo "⚠️ $(TIMESTAMP) Tải SnpEff DB ${SNPEFF_DB}..."
  java ${JAVA_OPTS_SNPEFF} -jar "$SNPEFF_JAR" download "$SNPEFF_DB" -c "$SNPEFF_CONFIG"
fi

# [1] SnpEff
java ${JAVA_OPTS_SNPEFF} -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
  -stats "${SNPEFF_STATS}" \
  "${DV_VCF}" | bgzip -@ ${THREADS} -c > "${ANN_VCF_SNPEFF}"
tabix -f -p vcf "${ANN_VCF_SNPEFF}"

# [2] gnomAD
[ -f "${GNOMAD_VCF}" ] || { echo "❌ Không thấy gnomAD: ${GNOMAD_VCF}"; exit 1; }
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -id -info AF \
  "${GNOMAD_VCF}" "${ANN_VCF_SNPEFF}" \
  | bgzip -@ ${THREADS} -c > "${ANN_VCF_GNOMAD}"
tabix -f -p vcf "${ANN_VCF_GNOMAD}"

# Đổi AF → GNOMAD_AF (nếu có bcftools)
if command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
  echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="Allele frequency from gnomAD">' > "${ANN_DIR}/gn_hdr.hdr"
  ${BCFTOOLS_BIN} annotate -h "${ANN_DIR}/gn_hdr.hdr" \
    -c INFO/GNOMAD_AF:=INFO/AF -x INFO/AF \
    -O z -o "${ANN_DIR}/_tmp_gnomad_renamed.vcf.gz" "${ANN_VCF_GNOMAD}"
  tabix -f -p vcf "${ANN_DIR}/_tmp_gnomad_renamed.vcf.gz"
  mv -f "${ANN_DIR}/_tmp_gnomad_renamed.vcf.gz" "${ANN_VCF_GNOMAD}"
  mv -f "${ANN_DIR}/_tmp_gnomad_renamed.vcf.gz.tbi" "${ANN_VCF_GNOMAD}.tbi"
  rm -f "${ANN_DIR}/gn_hdr.hdr"
fi

# [3] ClinVar
[ -f "${CLINVAR_VCF}" ] || { echo "❌ Không thấy ClinVar: ${CLINVAR_VCF}"; exit 1; }
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info CLNSIG,CLNDN \
  "${CLINVAR_VCF}" "${ANN_VCF_GNOMAD}" \
  | bgzip -@ ${THREADS} -c > "${ANN_VCF_CLINVAR}"
tabix -f -p vcf "${ANN_VCF_CLINVAR}"

# [4] 1000 Genomes (+ đổi AF → KG_AF nếu có bcftools)
[ -f "${THOUSANDG_VCF}" ] || { echo "❌ Không thấy 1000G: ${THOUSANDG_VCF}"; exit 1; }
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info AF \
  "${THOUSANDG_VCF}" "${ANN_VCF_CLINVAR}" > "${ANN_DIR}/_tmp_kg.vcf"

if command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1; then
  echo '##INFO=<ID=KG_AF,Number=A,Type=Float,Description="Allele frequency from 1000 Genomes Project">' > "${ANN_DIR}/kg_hdr.hdr"
  ${BCFTOOLS_BIN} annotate -h "${ANN_DIR}/kg_hdr.hdr" \
    -c INFO/KG_AF:=INFO/AF -x INFO/AF \
    -O z -o "${ANN_VCF_FINAL}" "${ANN_DIR}/_tmp_kg.vcf"
  tabix -f -p vcf "${ANN_VCF_FINAL}"
  rm -f "${ANN_DIR}/kg_hdr.hdr" "${ANN_DIR}/_tmp_kg.vcf"
else
  bgzip -f -@ ${THREADS} -c "${ANN_DIR}/_tmp_kg.vcf" > "${ANN_VCF_FINAL}"
  tabix -f -p vcf "${ANN_VCF_FINAL}"
  rm -f "${ANN_DIR}/_tmp_kg.vcf"
fi

echo "✅ Annotated VCF cuối: ${ANN_VCF_FINAL}"

# ===================================================================
# 5) COVERAGE
# ===================================================================
echo "---=== BƯỚC 5: Coverage (mosdepth) ===---"
conda activate "${ENV_BRCA}"
${MOSDEPTH_BIN} --threads ${THREADS} -n --by ${TARGET_BED} "${COVERAGE_PREFIX}" "${RECAL_BAM}"

# ===================================================================
# 6) MULTIQC
# ===================================================================
echo "---=== BƯỚC 6: MultiQC ===---"
conda activate "${ENV_MQC}"
${MULTIQC_BIN} "${SAMPLE_DIR}" \
  --outdir "${MULTIQC_DIR}" \
  --title "Báo cáo QC cho mẫu ${SAMPLE_NAME}" \
  --filename "${MULTIQC_FILENAME}" \
  --force

# ===================================================================
# 7) DỌN DẸP
# ===================================================================
if [ "$CLEANUP" = true ]; then
  echo "---=== BƯỚC 7: Dọn dẹp ===---"
  rm -f "${TRIMMED_R1}" "${UNPAIRED_R1}" \
        "${TRIMMED_R2}" "${UNPAIRED_R2}" \
        "${BAM}" "${SORTED_BAM}" \
        "${BAM_DEDUP}" "${BAM_DEDUP_BAI}" \
        "${RECAL_DATA_TABLE}"
fi

echo ""
echo "==================================================================="
echo "==  PIPELINE (DeepVariant) HOÀN TẤT CHO MẪU: ${SAMPLE_NAME}"
echo "==  Kết quả & báo cáo tại: ${SAMPLE_DIR}"
echo "==  DeepVariant VCF: ${DV_VCF}"
echo "==  Annotated VCF  : ${ANN_VCF_FINAL}"
echo "==  MultiQC        : ${MULTIQC_DIR}/${MULTIQC_FILENAME}"
echo "==================================================================="
