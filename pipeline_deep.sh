#!/usr/bin/env bash
# ===================================================================
# ==      PIPELINE PHÂN TÍCH DỮ LIỆU BRCA — DeepVariant (1 SM)    ==
# ==      Có cập nhật SnpEff/SnpSift theo khối annotate mới       ==
# ===================================================================
set -euo pipefail

# ======================= THAM SỐ & CẤU HÌNH =======================
SAMPLE_NAME="${1:?Vui lòng truyền SAMPLE_NAME}"
THREADS="${THREADS:-8}"
CLEANUP="${CLEANUP:-true}"
TIMESTAMP() { date '+%Y-%m-%d %H:%M:%S'; }

# Conda envs
ENV_BRCA="${ENV_BRCA:-BRCA}"     # fastqc, trimmomatic, bwa, samtools, mosdepth
ENV_GATK="${ENV_GATK:-GATK}"     # gatk (BQSR)
ENV_MQC="${ENV_MQC:-MQC}"        # multiqc
ENV_ANN="${ENV_ANN:-BCF}"        # bcftools, vt (optional), tabix/bgzip, java

# Project + log
PROJECT_DIR="${PROJECT_DIR:-/media/shmily/writable/BRCA_project}"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_DIR}"

# Tham chiếu & mục tiêu
REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Known sites cho BQSR
KNOWN_SNP="${PROJECT_DIR}/reference/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDEL="${PROJECT_DIR}/reference/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS_1000G_INDEL="${PROJECT_DIR}/reference/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# Tools
FASTQC_BIN="${FASTQC_BIN:-fastqc}"
TRIMMOMATIC_BIN="${TRIMMOMATIC_BIN:-trimmomatic}"
BWA_BIN="${BWA_BIN:-bwa}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
MOSDEPTH_BIN="${MOSDEPTH_BIN:-mosdepth}"
MULTIQC_BIN="${MULTIQC_BIN:-multiqc}"
TABIX_BIN="${TABIX_BIN:-tabix}"
BGZIP_BIN="${BGZIP_BIN:-bgzip}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"

# DeepVariant (Docker image)
DV_DOCKER_IMAGE="${DV_DOCKER_IMAGE:-google/deepvariant:1.9.0}"

# Picard + snpEff/SnpSift
PICARD_JAR="${PICARD_JAR:-/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar}"
SNPEFF_HOME="${PROJECT_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
SNPEFF_DB="${SNPEFF_DB:-GRCh38.86}"
SNPEFF_DATA_DIR="${SNPEFF_HOME}/data"

JAVA_OPTS_SNPEFF="-Xmx6g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_SNPSIFT="-Xmx4g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_PICARD="-Xmx4g -Djava.awt.headless=true -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"

# External annotation
GNOMAD_VCF="${PROJECT_DIR}/reference/resources/gnomad.v4.1.panel.merged.vcf.gz"
CLINVAR_VCF="${PROJECT_DIR}/reference/resources/clinvar_20250810.vcf.gz"
THOUSANDG_VCF="${PROJECT_DIR}/reference/resources/1000g.panel.merged.vcf.gz"

# Input FASTQ
RAW_DIR="${PROJECT_DIR}/raw_data"
READ1="${RAW_DIR}/${SAMPLE_NAME}_1.fastq.gz"
READ2="${RAW_DIR}/${SAMPLE_NAME}_2.fastq.gz"
ADAPTER_FILE="${ADAPTER_FILE:-/home/shmily/miniconda/envs/BRCA/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa}"

# Thư mục bước
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
FASTQC_RAW_DIR="${SAMPLE_DIR}/fastqc_raw"
FASTQC_TRIM_DIR="${SAMPLE_DIR}/fastqc_trimmed"
BWA_DIR="${SAMPLE_DIR}/Bwa_alignments"
RECAL_DIR="${SAMPLE_DIR}/recal"
DEEP_DIR="${SAMPLE_DIR}/Deepvariants"
ANN_DIR="${SAMPLE_DIR}/snpeff"
COVERAGE_DIR="${SAMPLE_DIR}/coverage"
MULTIQC_DIR="${SAMPLE_DIR}/multiqc_report"
mkdir -p "${TRIM_DIR}" "${FASTQC_RAW_DIR}" "${FASTQC_TRIM_DIR}" \
         "${BWA_DIR}" "${RECAL_DIR}" "${DEEP_DIR}" "${ANN_DIR}" \
         "${COVERAGE_DIR}" "${MULTIQC_DIR}"

# File trung gian & đầu ra
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"
UNPAIRED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_unpaired.fastq.gz"
UNPAIRED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_unpaired.fastq.gz"

BAM="${BWA_DIR}/${SAMPLE_NAME}_aligned.bam"
SORTED_BAM="${BWA_DIR}/${SAMPLE_NAME}_sorted.bam"
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam}"
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"   # sửa lỗi ngoặc
BAM_DEDUP_BAI="${BWA_DIR}/${SAMPLE_NAME}_dedup.bai"
DEDUP_METRICS="${BWA_DIR}/${SAMPLE_NAME}_dedup_metrics.txt"
STATS_FILE="${BWA_DIR}/${SAMPLE_NAME}_samtools_stats.txt"
FLAGSTAT_FILE="${BWA_DIR}/${SAMPLE_NAME}_samtools_flagstat.txt"

RECAL_DATA_TABLE="${RECAL_DIR}/${SAMPLE_NAME}_recal_data.table"
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

DV_VCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.vcf.gz"
DV_GVCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.g.vcf.gz"

# Annotate outputs (đặt hậu tố _dv_)
ANN_SNPEFF_VCF="${ANN_DIR}/${SAMPLE_NAME}_dv_snpeff.vcf"
ANN_GNOMAD_VCF="${ANN_DIR}/${SAMPLE_NAME}_dv_gnomad.vcf"
ANN_GNOMAD_RENAMED="${ANN_DIR}/${SAMPLE_NAME}_dv_gnomad_renamed.vcf"
ANN_CLINVAR_VCF="${ANN_DIR}/${SAMPLE_NAME}_dv_clinvar.vcf"
ANN_1KG_VCF="${ANN_DIR}/${SAMPLE_NAME}_dv_1kg.vcf"
ANN_VCF_FINAL="${ANN_DIR}/${SAMPLE_NAME}_dv_final_annotated.vcf"

# Các tệp tạm trong khối annotate mới
NORM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.norm.vcf.gz"
ATOM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.atom.vcf.gz"
HARM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.harm.vcf.gz"

# Coverage & report
COVERAGE_PREFIX="${COVERAGE_DIR}/${SAMPLE_NAME}"
MULTIQC_FILENAME="${SAMPLE_NAME}_dv_report.html"

# Read group
RG_ID="${SAMPLE_NAME}_RG"; PLATFORM="Illumina"; LIBRARY_ID="Lib1"
PLATFORM_UNIT="${SAMPLE_NAME}_${PLATFORM}_${LIBRARY_ID}"
RG_STR="@RG\tID:${RG_ID}\tSM:${SAMPLE_NAME}\tPL:${PLATFORM}\tLB:${LIBRARY_ID}\tPU:${PLATFORM_UNIT}"

# ======================= HÀM HỖ TRỢ ANNOTATE ======================
need_cmd()  { command -v "$1" >/dev/null 2>&1 || { echo "❌ Thiếu tool: $1"; exit 1; }; }
need_file() { [[ -f "$1" ]] || { echo "❌ Thiếu file: $1"; exit 1; }; }

supports_bcftools_atomize() {
  local exe; exe="$(command -v "${BCFTOOLS_BIN}" || true)"; [[ -n "$exe" ]] || return 1
  local help; help="$("$exe" norm -h 2>&1 || true)"
  case "$help" in *"--atomize"*) return 0 ;; *) return 1 ;; esac
}

has_chr() {
  "${BCFTOOLS_BIN}" view -h "$1" | grep -m1 '^##contig' | grep -q 'ID=chr' && return 0 || return 1
}

MKMAP_ADDCHR="$(mktemp)"; cat > "$MKMAP_ADDCHR" <<'EOF'
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
MKMAP_RMCHR="$(mktemp)"; awk '{print $2"\t"$1}' "$MKMAP_ADDCHR" > "$MKMAP_RMCHR"

harmonize_to_resource() {
  local sample="$1" resource="$2" out="$3"
  if has_chr "$sample" && ! has_chr "$resource"; then
    echo "[$(TIMESTAMP)] 🔁 Bỏ 'chr' để khớp $(basename "$resource")"
    "${BCFTOOLS_BIN}" annotate --rename-chrs "$MKMAP_RMCHR" -O z -o "$out" "$sample"
    "${TABIX_BIN}" -f -p vcf "$out"
  elif ! has_chr "$sample" && has_chr "$resource"; then
    echo "[$(TIMESTAMP)] 🔁 Thêm 'chr' để khớp $(basename "$resource")"
    "${BCFTOOLS_BIN}" annotate --rename-chrs "$MKMAP_ADDCHR" -O z -o "$out" "$sample"
    "${TABIX_BIN}" -f -p vcf "$out"
  else
    cp -f "$sample" "$out"
    "${TABIX_BIN}" -f -p vcf "$out" || true
  fi
}

rename_AF_to_new_tag() {
  local in="$1" newtag="$2" out="$3"
  local hdr tmp tsv_gz
  hdr="$(mktemp)"
  echo "##INFO=<ID=${newtag},Number=A,Type=Float,Description=\"Allele frequency from ${newtag}\">" > "$hdr"
  tmp="$(mktemp)"
  "${BCFTOOLS_BIN}" query -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$in" > "$tmp" || true
  tsv_gz="${tmp}.gz"; "${BGZIP_BIN}" -f -c "$tmp" > "$tsv_gz" && "${TABIX_BIN}" -f -s 1 -b 2 -e 2 "$tsv_gz"
  "${BCFTOOLS_BIN}" annotate -a "$tsv_gz" -c CHROM,POS,REF,ALT,INFO/"$newtag" -h "$hdr" -x INFO/AF -O v -o "$out" "$in"
  rm -f "$hdr" "$tmp" "$tsv_gz" "${tsv_gz}.tbi"
}

ensure_ref_index() {
  if [[ ! -f "${REF}.fai" ]]; then
    echo "[$(TIMESTAMP)] 🧩 Tạo index FASTA tham chiếu (.fai)..."
    "${SAMTOOLS_BIN}" faidx "$REF"
  fi
}

# ======================= KÍCH HOẠT CONDA & CHECK ===================
eval "$(conda shell.bash hook)"
for f in "${READ1}" "${READ2}" "${ADAPTER_FILE}" "${REF}" "${TARGET_BED}"; do need_file "$f"; done

# ===================================================================
# 1) QC & TRIMMING
# ===================================================================
echo -e "\n---=== BƯỚC 1: QC & Trimming ===---"
conda activate "${ENV_BRCA}"
"${FASTQC_BIN}" --threads "${THREADS}" -o "${FASTQC_RAW_DIR}" "${READ1}" "${READ2}"
"${TRIMMOMATIC_BIN}" PE -threads "${THREADS}" -phred33 \
  "${READ1}" "${READ2}" \
  "${TRIMMED_R1}" "${UNPAIRED_R1}" \
  "${TRIMMED_R2}" "${UNPAIRED_R2}" \
  ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
"${FASTQC_BIN}" --threads "${THREADS}" -o "${FASTQC_TRIM_DIR}" "${TRIMMED_R1}" "${TRIMMED_R2}"

# ===================================================================
# 2) ALIGN + DEDUP
# ===================================================================
echo -e "\n---=== BƯỚC 2: Gióng hàng & đánh dấu trùng lặp ===---"
"${BWA_BIN}" mem -Y -K 100000000 -t "${THREADS}" -R "${RG_STR}" "${REF}" "${TRIMMED_R1}" "${TRIMMED_R2}" | \
  "${SAMTOOLS_BIN}" view -b -o "${BAM}" -
java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" SortSam I="${BAM}" O="${SORTED_BAM}" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" MarkDuplicates I="${SORTED_BAM}" O="${BAM_DEDUP}" M="${DEDUP_METRICS}" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
"${SAMTOOLS_BIN}" stats "${BAM_DEDUP}" > "${STATS_FILE}"
"${SAMTOOLS_BIN}" flagstat "${BAM_DEDUP}" > "${FLAGSTAT_FILE}"

# ===================================================================
# 3) BQSR
# ===================================================================
echo -e "\n---=== BƯỚC 3: BQSR ===---"
conda activate "${ENV_GATK}"
ensure_ref_index
gatk BaseRecalibrator -I "${BAM_DEDUP}" -R "${REF}" --known-sites "${KNOWN_SNP}" --known-sites "${KNOWN_INDEL}" --known-sites "${MILLS_1000G_INDEL}" -O "${RECAL_DATA_TABLE}"
gatk ApplyBQSR -R "${REF}" -I "${BAM_DEDUP}" --bqsr-recal-file "${RECAL_DATA_TABLE}" -O "${RECAL_BAM}"
"${SAMTOOLS_BIN}" index "${RECAL_BAM}"

# ===================================================================
# 4) DEEPVARIANT + ANNOTATE (MỚI)
# ===================================================================
echo -e "\n---=== BƯỚC 4: DeepVariant & Annotate (mới) ===---"
docker run --rm -u "$(id -u)":"$(id -g)" -v "${PROJECT_DIR}:${PROJECT_DIR}" "${DV_DOCKER_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref="${REF}" \
    --reads="${RECAL_BAM}" \
    --regions="${TARGET_BED}" \
    --output_vcf="${DV_VCF}" \
    --output_gvcf="${DV_GVCF}" \
    --num_shards="${THREADS}" \
    --vcf_stats_report=true
"${TABIX_BIN}" -f -p vcf "${DV_VCF}" || true
"${TABIX_BIN}" -f -p vcf "${DV_GVCF}" || true

# ---- Khối annotate mới ----
conda activate "${ENV_ANN}"
for c in "${BCFTOOLS_BIN}" "${TABIX_BIN}" "${BGZIP_BIN}" java; do need_cmd "$c"; done
need_file "${GNOMAD_VCF}"; need_file "${CLINVAR_VCF}"; need_file "${THOUSANDG_VCF}"
[[ -d "${SNPEFF_DATA_DIR}/${SNPEFF_DB}" ]] || java -jar "${SNPEFF_JAR}" download "${SNPEFF_DB}" -c "${SNPEFF_CONFIG}"

# Normalize + split multi-allelic
echo "[$(TIMESTAMP)] 🔧 Normalize + split..."
"${BCFTOOLS_BIN}" norm -m -both -f "${REF}" -O z -o "${NORM_VCF_GZ}" "${DV_VCF}"
"${TABIX_BIN}" -f -p vcf "${NORM_VCF_GZ}"

# Atomize
ATOMIZE="${ATOMIZE:-true}"
if [[ "${ATOMIZE}" == "true" ]] && supports_bcftools_atomize; then
  echo "[$(TIMESTAMP)] 🧩 Atomize bằng bcftools --atomize..."
  "${BCFTOOLS_BIN}" norm --atomize -f "${REF}" -O z -o "${ATOM_VCF_GZ}" "${NORM_VCF_GZ}"
  "${TABIX_BIN}" -f -p vcf "${ATOM_VCF_GZ}"
elif [[ "${ATOMIZE}" == "true" ]] && command -v vt >/dev/null 2>&1; then
  echo "[$(TIMESTAMP)] 🧩 Atomize bằng vt decompose -s..."
  tmp="${ATOM_VCF_GZ%.gz}"
  "${BCFTOOLS_BIN}" view -Ov -o "$tmp" "${NORM_VCF_GZ}"
  vt decompose -s "$tmp" -o "${tmp%.vcf}.atom.vcf"
  "${BGZIP_BIN}" -f "${tmp%.vcf}.atom.vcf"
  "${TABIX_BIN}" -f -p vcf "${tmp%.vcf}.atom.vcf.gz"
  mv -f "${tmp%.vcf}.atom.vcf.gz" "${ATOM_VCF_GZ}"
  mv -f "${tmp%.vcf}.atom.vcf.gz.tbi" "${ATOM_VCF_GZ}.tbi"
  rm -f "$tmp"
else
  echo "[$(TIMESTAMP)] ⏭️  Bỏ qua atomize."
  cp -f "${NORM_VCF_GZ}" "${ATOM_VCF_GZ}"; "${TABIX_BIN}" -f -p vcf "${ATOM_VCF_GZ}" || true
fi

# Harmonize chr với gnomAD
echo "[$(TIMESTAMP)] 🔁 Harmonize chr với gnomAD..."
harmonize_to_resource "${ATOM_VCF_GZ}" "${GNOMAD_VCF}" "${HARM_VCF_GZ}"

# SnpEff
echo "[$(TIMESTAMP)] 🔬 SnpEff..."
java ${JAVA_OPTS_SNPEFF} -jar "${SNPEFF_JAR}" ann -c "${SNPEFF_CONFIG}" -v "${SNPEFF_DB}" -stats "${ANN_DIR}/${SAMPLE_NAME}_dv" "${HARM_VCF_GZ}" > "${ANN_SNPEFF_VCF}"

# gnomAD (AF) → GNOMAD_AF
echo "[$(TIMESTAMP)] 📊 gnomAD annotate..."
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info AF "${GNOMAD_VCF}" "${ANN_SNPEFF_VCF}" > "${ANN_GNOMAD_VCF}"
echo "[$(TIMESTAMP)] 📝 Đổi AF → GNOMAD_AF..."
rename_AF_to_new_tag "${ANN_GNOMAD_VCF}" "GNOMAD_AF" "${ANN_GNOMAD_RENAMED}"

# ClinVar
echo "[$(TIMESTAMP)] 🧬 ClinVar annotate..."
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info CLNSIG,CLNDN "${CLINVAR_VCF}" "${ANN_GNOMAD_RENAMED}" > "${ANN_CLINVAR_VCF}"

# 1000 Genomes (AF) → KG_AF
echo "[$(TIMESTAMP)] 🌍 1000G annotate..."
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info AF "${THOUSANDG_VCF}" "${ANN_CLINVAR_VCF}" > "${ANN_1KG_VCF}"
echo "[$(TIMESTAMP)] 📝 Đổi AF → KG_AF (final)..."
rename_AF_to_new_tag "${ANN_1KG_VCF}" "KG_AF" "${ANN_VCF_FINAL}"

# (Tuỳ chọn) nén cuối + index
if command -v "${BGZIP_BIN}" >/dev/null 2>&1; then
  "${BGZIP_BIN}" -f "${ANN_VCF_FINAL}"
  "${TABIX_BIN}" -f -p vcf "${ANN_VCF_FINAL}.gz"
fi

# ===================================================================
# 5) COVERAGE
# ===================================================================
echo -e "\n---=== BƯỚC 5: Mosdepth ===---"
conda activate "${ENV_BRCA}"
"${MOSDEPTH_BIN}" --threads "${THREADS}" -n --by "${TARGET_BED}" "${COVERAGE_PREFIX}" "${RECAL_BAM}"

# ===================================================================
# 6) MULTIQC
# ===================================================================
echo -e "\n---=== BƯỚC 6: MultiQC ===---"
conda activate "${ENV_MQC}"
"${MULTIQC_BIN}" "${SAMPLE_DIR}" --outdir "${MULTIQC_DIR}" --title "Báo cáo QC cho mẫu ${SAMPLE_NAME}" --filename "${MULTIQC_FILENAME}" --force

# ===================================================================
# 7) CLEANUP
# ===================================================================
echo -e "\n---=== BƯỚC 7: Dọn dẹp ===---"
if [[ "${CLEANUP}" == "true" ]]; then
  rm -f \
    "${TRIMMED_R1}" "${UNPAIRED_R1}" "${TRIMMED_R2}" "${UNPAIRED_R2}" \
    "${BAM}" "${SORTED_BAM}" "${BAM_DEDUP}" "${BAM_DEDUP_BAI}" \
    "${RECAL_DATA_TABLE}" \
    "${DV_GVCF}" "${DV_GVCF}.tbi" \
    "${DV_VCF}.stats.html" 2>/dev/null || true \
    "${NORM_VCF_GZ}" "${NORM_VCF_GZ}.tbi" \
    "${ATOM_VCF_GZ}" "${ATOM_VCF_GZ}.tbi" \
    "${HARM_VCF_GZ}" "${HARM_VCF_GZ}.tbi" \
    "${ANN_SNPEFF_VCF}" \
    "${ANN_GNOMAD_VCF}" "${ANN_GNOMAD_RENAMED}" \
    "${ANN_CLINVAR_VCF}" "${ANN_1KG_VCF}" \
    "$MKMAP_ADDCHR" "$MKMAP_RMCHR"
else
  echo "Bỏ qua dọn dẹp (CLEANUP=false)"
fi

echo -e "\n==================================================================="
echo "✅ HOÀN TẤT (DeepVariant) — ${SAMPLE_NAME}"
echo "Annotated VCF: ${ANN_VCF_FINAL} (và .gz nếu đã nén)"
echo "Báo cáo MultiQC: ${MULTIQC_DIR}/${MULTIQC_FILENAME}"
echo "==================================================================="
