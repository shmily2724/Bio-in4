#!/bin/bash
# ===================================================================
# == BRCA PIPELINE (DeepVariant) — TOOL-CHECK ONLY — 2025-08-16    ==
# ===================================================================
set -euo pipefail

# ========== CẤU HÌNH CHUNG ==========
SAMPLE_NAME="${1:?Vui lòng truyền SAMPLE_NAME}"
THREADS="${THREADS:-8}"
CLEANUP="${CLEANUP:-true}"
TIMESTAMP() { date '+%Y-%m-%d %H:%M:%S'; }

ENV_BRCA="BRCA"
ENV_GATK="GATK"
ENV_MQC="MQC"

PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
LOG="${SAMPLE_DIR}/${SAMPLE_NAME}_pipeline.log"
mkdir -p "${SAMPLE_DIR}"
# exec > >(tee -i "$LOG") 2>&1

# ========== THAM CHIẾU & VÙNG MỤC TIÊU ==========
REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Known sites cho BQSR (KHÔNG kiểm tra file/index)
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

PICARD_JAR="/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar"

# ========== SnpEff/SnpSift (KHÔNG kiểm tra DB) ==========
SNPEFF_HOME="${PROJECT_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
SNPEFF_DB="GRCh38.86"

JAVA_OPTS_SNPEFF="-Xmx6g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_SNPSIFT="-Xmx4g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_PICARD="-Xmx4g -Djava.awt.headless=true -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"

# ========== EXTERNAL RESOURCES (KHÔNG kiểm tra) ==========
GNOMAD_VCF="${PROJECT_DIR}/reference/resources/gnomad.v4.1.panel.merged.vcf.gz"
CLINVAR_VCF="${PROJECT_DIR}/reference/resources/clinvar_20250810.vcf.gz"
THOUSANDG_VCF="${PROJECT_DIR}/reference/resources/1000g.panel.merged.vcf.gz"

# ========== THƯ MỤC KẾT QUẢ ==========
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
FASTQC_RAW_DIR="${SAMPLE_DIR}/fastqc_raw"
FASTQC_TRIM_DIR="${SAMPLE_DIR}/fastqc_trimmed"
BWA_DIR="${SAMPLE_DIR}/Bwa_alignments"
RECAL_DIR="${SAMPLE_DIR}/recal"
DEEP_DIR="${SAMPLE_DIR}/Deepvariants"
ANN_DIR="${SAMPLE_DIR}/snpeff"
COVERAGE_DIR="${SAMPLE_DIR}/coverage"
MULTIQC_DIR="${SAMPLE_DIR}/multiqc_report"
mkdir -p "${TRIM_DIR}" "${FASTQC_RAW_DIR}" "${FASTQC_TRIM_DIR}" "${BWA_DIR}" \
         "${RECAL_DIR}" "${DEEP_DIR}" "${ANN_DIR}" "${COVERAGE_DIR}" "${MULTIQC_DIR}"

# ========== FILE TRUNG GIAN & ĐẦU RA ==========
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

# DeepVariant outputs
DV_VCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.vcf.gz"
DV_GVCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.g.vcf.gz"
DV_DOCKER_IMAGE="${DV_DOCKER_IMAGE:-google/deepvariant:1.9.0}"
DV_MODEL="${DV_MODEL:-WES}"   # WES|WGS|PACBIO|HYBRID_PACBIO_ILLUMINA

# Annotation outputs
SNPEFF_STATS_HTML="${ANN_DIR}/${SAMPLE_NAME}_dv_snpeff_stats.html"
ANN_VCF_FINAL="${ANN_DIR}/${SAMPLE_NAME}_dv_final_annotated.vcf"
ANN_VCF_SNPEFF="${ANN_DIR}/${SAMPLE_NAME}_dv_snpeff.vcf"
ANN_VCF_GNOMAD="${ANN_DIR}/${SAMPLE_NAME}_dv_gnomad.vcf"
ANN_VCF_CLINVAR="${ANN_DIR}/${SAMPLE_NAME}_dv_clinvar.vcf"

TMP1="${ANN_VCF_SNPEFF}"
TMP2="${ANN_VCF_GNOMAD}"
TMP3="${ANN_VCF_CLINVAR}"
ANN_VCF="${ANN_VCF_FINAL}"

COVERAGE_PREFIX="${COVERAGE_DIR}/${SAMPLE_NAME}"
MULTIQC_FILENAME="${SAMPLE_NAME}_dv_report.html"

RG_ID="${SAMPLE_NAME}_RG"
PLATFORM="Illumina"
LIBRARY_ID="Lib1"
PLATFORM_UNIT="${SAMPLE_NAME}_${PLATFORM}_${LIBRARY_ID}"
RG_STR="@RG\tID:${RG_ID}\tSM:${SAMPLE_NAME}\tPL:${PLATFORM}\tLB:${LIBRARY_ID}\tPU:${PLATFORM_UNIT}"

# --- KÍCH HOẠT CONDA ---
eval "$(conda shell.bash hook)"

# --- CHỈ KIỂM TRA TOOL (không kiểm tra database/file dữ liệu) ---
need() { command -v "$1" >/dev/null 2>&1 || { echo "❌ Thiếu tool: $1"; exit 127; }; }
need fastqc; need trimmomatic; need bwa; need samtools; need gatk; need mosdepth; need multiqc
need java; need bcftools; need bgzip; need tabix
# DeepVariant: cần Docker hoặc binary native
if command -v docker >/dev/null 2>&1; then
  USE_DOCKER=1
else
  need run_deepvariant
  USE_DOCKER=0
fi

# ===================================================================
# BƯỚC 1: QC & TRIMMING
# ===================================================================
printf '\n'; echo "---=== BƯỚC 1: QC & Trimming ===---"
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
# BƯỚC 2: ALIGN → SORT → MARKDUP
# ===================================================================
printf '\n'; echo "---=== BƯỚC 2: Align/Sort/MarkDup ===---"
${BWA_BIN} mem -Y -K 100000000 -t ${THREADS} -R "${RG_STR}" "$REF" "${TRIMMED_R1}" "${TRIMMED_R2}" \
  | ${SAMTOOLS_BIN} view -b -o "${BAM}" -
java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" SortSam I="${BAM}" O="${SORTED_BAM}" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" MarkDuplicates I="${SORTED_BAM}" O="${BAM_DEDUP}" M="${DEDUP_METRICS}" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
${SAMTOOLS_BIN} stats "${BAM_DEDUP}" > "${STATS_FILE}"; ${SAMTOOLS_BIN} flagstat "${BAM_DEDUP}" > "${FLAGSTAT_FILE}"

# ===================================================================
# BƯỚC 3: BQSR
# ===================================================================
printf '\n'; echo "---=== BƯỚC 3: BQSR ===---"
conda activate "${ENV_GATK}"
${GATK_BIN} BaseRecalibrator -I "${BAM_DEDUP}" -R "${REF}" \
  --known-sites "${KNOWN_SNP}" --known-sites "${KNOWN_INDEL}" --known-sites "${MILLS_1000G_INDEL}" \
  -O "${RECAL_DATA_TABLE}"
${GATK_BIN} ApplyBQSR -R "${REF}" -I "${BAM_DEDUP}" --bqsr-recal-file "${RECAL_DATA_TABLE}" -O "${RECAL_BAM}"
${SAMTOOLS_BIN} index "${RECAL_BAM}"

# ===================================================================
# BƯỚC 4: DEEPVARIANT → ANNOTATE
# ===================================================================
printf '\n'; echo "---=== BƯỚC 4: DeepVariant + Chú giải ===---"

if [[ "${USE_DOCKER}" -eq 1 ]]; then
  docker run --rm -v "${PROJECT_DIR}:${PROJECT_DIR}" -v "/etc/localtime:/etc/localtime:ro" \
    "${DV_DOCKER_IMAGE}" \
    /opt/deepvariant/bin/run_deepvariant \
      --model_type="${DV_MODEL}" \
      --ref="${REF}" \
      --reads="${RECAL_BAM}" \
      --regions="${TARGET_BED}" \
      --output_vcf="${DV_VCF}" \
      --output_gvcf="${DV_GVCF}" \
      --num_shards="${THREADS}"
else
  run_deepvariant --model_type="${DV_MODEL}" --ref="${REF}" --reads="${RECAL_BAM}" \
    --regions="${TARGET_BED}" --output_vcf="${DV_VCF}" --output_gvcf="${DV_GVCF}" --num_shards="${THREADS}"
fi

(tabix -f -p vcf "${DV_VCF}"  || true)
(tabix -f -p vcf "${DV_GVCF}" || true)

supports_atomize() { bcftools norm -h 2>&1 | grep -q -- '--atomize'; }

NORM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.norm.vcf.gz"
bcftools norm -m -both -f "${REF}" -O z -o "${NORM_VCF_GZ}" "${DV_VCF}"
tabix -f -p vcf "${NORM_VCF_GZ}"

ATOM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.atom.vcf.gz"
if supports_atomize; then
  bcftools norm --atomize -f "${REF}" -O z -o "${ATOM_VCF_GZ}" "${NORM_VCF_GZ}"
  tabix -f -p vcf "${ATOM_VCF_GZ}"
elif command -v vt >/dev/null 2>&1; then
  bcftools view -Ov -o "${ANN_DIR}/tmp.dv.norm.vcf" "${NORM_VCF_GZ}"
  vt decompose -s "${ANN_DIR}/tmp.dv.norm.vcf" -o "${ANN_DIR}/tmp.dv.atom.vcf"
  bgzip -f "${ANN_DIR}/tmp.dv.atom.vcf"; tabix -f -p vcf "${ANN_DIR}/tmp.dv.atom.vcf.gz"
  mv -f "${ANN_DIR}/tmp.dv.atom.vcf.gz" "${ATOM_VCF_GZ}"
  mv -f "${ANN_DIR}/tmp.dv.atom.vcf.gz.tbi" "${ATOM_VCF_GZ}.tbi"
  rm -f "${ANN_DIR}/tmp.dv.norm.vcf"
else
  echo "⚠️ Không có --atomize hoặc vt → bỏ qua atomize."
  cp -f "${NORM_VCF_GZ}" "${ATOM_VCF_GZ}"; tabix -f -p vcf "${ATOM_VCF_GZ}" || true
fi

ADDCHR_MAP="$(mktemp)"; cat > "$ADDCHR_MAP" <<'EOF'
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
X chrX
Y chrY
MT chrM
EOF
RMCHR_MAP="$(mktemp)"; awk '{print $2"\t"$1}' "$ADDCHR_MAP" > "$RMCHR_MAP"

HARM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.harm.vcf.gz"
if bcftools view -h "${ATOM_VCF_GZ}" | grep -m1 '^##contig' | grep -q 'ID=chr'; then
  if ! bcftools view -h "${GNOMAD_VCF}" | grep -m1 '^##contig' | grep -q 'ID=chr'; then
    bcftools annotate --rename-chrs "$RMCHR_MAP" -O z -o "${HARM_VCF_GZ}" "${ATOM_VCF_GZ}"
  else
    cp -f "${ATOM_VCF_GZ}" "${HARM_VCF_GZ}" && tabix -f -p vcf "${HARM_VCF_GZ}" || true
  fi
else
  if bcftools view -h "${GNOMAD_VCF}" | grep -m1 '^##contig' | grep -q 'ID=chr'; then
    bcftools annotate --rename-chrs "$ADDCHR_MAP" -O z -o "${HARM_VCF_GZ}" "${ATOM_VCF_GZ}"
  else
    cp -f "${ATOM_VCF_GZ}" "${HARM_VCF_GZ}" && tabix -f -p vcf "${HARM_VCF_GZ}" || true
  fi
fi
(tabix -f -p vcf "${HARM_VCF_GZ}" || true)
rm -f "$ADDCHR_MAP" "$RMCHR_MAP"

# SnpEff + SnpSift (KHÔNG kiểm tra DB)
java ${JAVA_OPTS_SNPEFF} -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
  -stats "${SNPEFF_STATS_HTML}" \
  "${HARM_VCF_GZ}" > "${TMP1}"

echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="Allele frequency from gnomAD">' > "${ANN_DIR}/gn_hdr.hdr"
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info AF "${GNOMAD_VCF}" "${TMP1}" > "${TMP2}"
bcftools annotate -h "${ANN_DIR}/gn_hdr.hdr" -c INFO/GNOMAD_AF:=INFO/AF -x INFO/AF -O v -o "${TMP2}.renamed" "${TMP2}"
mv -f "${TMP2}.renamed" "${TMP2}"

java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info CLNSIG,CLNDN "${CLINVAR_VCF}" "${TMP2}" > "${TMP3}"

echo '##INFO=<ID=KG_AF,Number=A,Type=Float,Description="Allele frequency from 1000 Genomes Project">' > "${ANN_DIR}/kg_hdr.hdr"
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info AF "${THOUSANDG_VCF}" "${TMP3}" > "${ANN_VCF}"
bcftools annotate -h "${ANN_DIR}/kg_hdr.hdr" -c INFO/KG_AF:=INFO/AF -x INFO/AF -O v -o "${ANN_VCF}.tmp" "${ANN_VCF}"
mv -f "${ANN_VCF}.tmp" "${ANN_VCF}"

bgzip -f "${ANN_VCF}"; tabix -f -p vcf "${ANN_VCF}.gz"
echo "✅ Annotated VCF (DV): ${ANN_VCF}.gz"

# ===================================================================
# BƯỚC 5: COVERAGE
# ===================================================================
printf '\n'; echo "---=== BƯỚC 5: Mosdepth ===---"
conda activate "${ENV_BRCA}"
${MOSDEPTH_BIN} --threads ${THREADS} -n --by ${TARGET_BED} "${COVERAGE_PREFIX}" "${RECAL_BAM}"

# ===================================================================
# BƯỚC 6: MULTIQC
# ===================================================================
printf '\n'; echo "---=== BƯỚC 6: MultiQC ===---"
conda activate "${ENV_MQC}"
${MULTIQC_BIN} "${SAMPLE_DIR}" --outdir "${MULTIQC_DIR}" --title "Báo cáo QC cho mẫu ${SAMPLE_NAME}" --filename "${MULTIQC_FILENAME}" --force

# ===================================================================
# BƯỚC 7: CLEANUP
# ===================================================================
if [ "${CLEANUP}" = true ]; then
  printf '\n'; echo "---=== BƯỚC 7: Cleanup ===---"
  rm -f \
    "${TRIMMED_R1}" "${UNPAIRED_R1}" "${TRIMMED_R2}" "${UNPAIRED_R2}" \
    "${BAM}" "${SORTED_BAM}" "${BAM_DEDUP}" "${BAM_DEDUP_BAI}" \
    "${RECAL_DATA_TABLE}" \
    "${DV_GVCF}" "${DV_GVCF}.tbi" \
    "${ANN_DIR}/tmp.dv.norm.vcf" "${ANN_DIR}/tmp.dv.atom.vcf" \
    "${ANN_DIR}/tmp.dv.atom.vcf.gz" "${ANN_DIR}/tmp.dv.atom.vcf.gz.tbi" \
    "${ANN_DIR}/${SAMPLE_NAME}_dv.norm.vcf.gz" "${ANN_DIR}/${SAMPLE_NAME}_dv.norm.vcf.gz.tbi" \
    "${ANN_DIR}/${SAMPLE_NAME}_dv.atom.vcf.gz" "${ANN_DIR}/${SAMPLE_NAME}_dv.atom.vcf.gz.tbi" \
    "${ANN_DIR}/${SAMPLE_NAME}_dv.harm.vcf.gz" "${ANN_DIR}/${SAMPLE_NAME}_dv.harm.vcf.gz.tbi" \
    "${TMP2}.renamed" "${ANN_VCF}.tmp" \
    "${ANN_DIR}/gn_hdr.hdr" "${ANN_DIR}/kg_hdr.hdr" \
    "${TMP1}" "${TMP2}" "${TMP3}" || true
fi

printf '\n'; echo "==================================================================="
echo "==  HOÀN TẤT DEEPVARIANT CHO MẪU: ${SAMPLE_NAME}"
echo "==  Kết quả: ${SAMPLE_DIR}"
echo "==  DeepVariant VCF: ${DV_VCF}"
echo "==  Annotated VCF : ${ANN_VCF}.gz"
echo "==  Báo cáo MultiQC: ${MULTIQC_DIR}/${MULTIQC_FILENAME}"
echo "==================================================================="
