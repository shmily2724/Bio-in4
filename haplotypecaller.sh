#!/bin/bash

# Tự động dừng script khi có lỗi
set -e

# Kiểm tra tham số đầu vào
if [ -z "$1" ]; then
    echo "Lỗi: Vui lòng cung cấp tên mẫu (SAMPLE_NAME)."
    echo "Cách dùng: $0 <sample_name>"
    exit 1
fi

SAMPLE_NAME=$1

# --- TẠO CÁC BIẾN ĐƯỜNG DẪN ---
PROJECT_DIR="/media/shmily/writable/BRCA_project"
REF="/media/shmily/writable/BRCA_project/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed" #link tải ở đây https://support.illumina.com/downloads/enrichment-bed-files-hg38.html
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
RECAL_DIR="${SAMPLE_DIR}/recal"
HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"

# Tạo các thư mục output
mkdir -p "${HAPLO_DIR}"

# File input
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

# File output
HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"

# Kích hoạt môi trường Conda (chỉ một lần)
eval "$(conda shell.bash hook)"
conda activate GATK

# --- BƯỚC 1: GỌI BIẾN THỂ VỚI HaplotypeCaller ---
echo "Bắt đầu HaplotypeCaller cho mẫu: ${SAMPLE_NAME}"
# Nếu chỉ phân tích 1 mẫu, có thể bỏ -ERC GVCF để ra thẳng VCF.
gatk HaplotypeCaller \
    -R "${REF}" \
    -I "${RECAL_BAM}" \
    -O "${HAPLO_GVCF}" \
    -L "${TARGET_BED}" \
    -ERC GVCF

echo "HaplotypeCaller hoàn tất."

# --- BƯỚC 2: CHUYỂN GVCF SANG VCF VỚI GenotypeGVCFs ---
echo "Bắt đầu GenotypeGVCFs cho mẫu: ${SAMPLE_NAME}"

gatk GenotypeGVCFs \
    -R "${REF}" \
    -V "${HAPLO_GVCF}" \
    -O "${HAPLO_VCF}" \
    -L "${TARGET_BED}"

echo "GenotypeGVCFs hoàn tất."

# --- BƯỚC 3: ĐÁNH INDEX CHO FILE VCF ---
echo "Đánh index cho file VCF..."
tabix -f -p vcf "${HAPLO_VCF}" # -f để overwrite tránh dừng script tại đây khi chạy lại
echo "Index hoàn tất cho mẫu: ${SAMPLE_NAME}!"
#!/bin/bash
# ===================================================================
# ==   HaplotypeCaller + GenotypeGVCFs + Annotate (SnpEff/SnpSift) ==
# ===================================================================
#fix cho code cũ 
#set -euo pipefail

# --- THAM SỐ ---
#SAMPLE_NAME="${1:?Lỗi: Vui lòng cung cấp tên mẫu (SAMPLE_NAME). Cách dùng: $0 <sample_name>}"
#CLEAN_INTERMEDIATE="${CLEAN_INTERMEDIATE:-true}"   # true|false: xóa file trung gian annotate

# --- ĐƯỜNG DẪN CỐ ĐỊNH ---
#PROJECT_DIR="/media/shmily/writable/BRCA_project"
#REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
#TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

#SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
#RECAL_DIR="${SAMPLE_DIR}/recal"
#HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"
#ANN_DIR="${SAMPLE_DIR}/snpeff"

#mkdir -p "${HAPLO_DIR}" "${ANN_DIR}"

# --- INPUT/OUTPUT CHO GỌI BIẾN THỂ ---
#RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"
#HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
#HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"

# --- RESOURCES CHO ANNOTATE ---
#GNOMAD_VCF="${PROJECT_DIR}/reference/resources/gnomad.v4.1.panel.merged.vcf.gz"
#CLINVAR_VCF="${PROJECT_DIR}/reference/resources/clinvar_20250810.vcf.gz"
#THOUSANDG_VCF="${PROJECT_DIR}/reference/resources/1000g.panel.merged.vcf.gz"

# --- SnpEff/SnpSift ---
#SNPEFF_HOME="${PROJECT_DIR}/snpEff"
#SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
#SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
#SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
#SNPEFF_DB="GRCh38.86"
#JAVA_OPTS_SNPEFF="-Xmx6g"
#JAVA_OPTS_SNPSIFT="-Xmx4g"

# --- OUTPUT ANNOTATE (hậu tố _hc_) ---
#SNPEFF_STATS="${ANN_DIR}/${SAMPLE_NAME}_hc"
#OUT_SNPEFF="${ANN_DIR}/${SAMPLE_NAME}_hc_snpeff.vcf"
#OUT_GNOMAD="${ANN_DIR}/${SAMPLE_NAME}_hc_gnomad.vcf"
#OUT_CLINVAR="${ANN_DIR}/${SAMPLE_NAME}_hc_clinvar.vcf"
#OUT_FINAL="${ANN_DIR}/${SAMPLE_NAME}_hc_final_annotated.vcf"

# --- KÍCH HOẠT CONDA (GATK) ---
#eval "$(conda shell.bash hook)"
#conda activate GATK

# =========================
# BƯỚC 1: HaplotypeCaller
# =========================
#echo "🎯 HaplotypeCaller → ${HAPLO_GVCF}"
#gatk HaplotypeCaller \
 # -R "${REF}" \
  #-I "${RECAL_BAM}" \
 # -O "${HAPLO_GVCF}" \
# -L "${TARGET_BED}" \
 # -ERC GVCF

# =========================
# BƯỚC 2: GenotypeGVCFs
# =========================
#echo "🧮 GenotypeGVCFs → ${HAPLO_VCF}"
#gatk GenotypeGVCFs \
#-R "${REF}" \
 # -V "${HAPLO_GVCF}" \
 # -O "${HAPLO_VCF}" \
  #-L "${TARGET_BED}"

#echo "🔖 Index VCF..."
#tabix -f -p vcf "${HAPLO_VCF}" || true

# =========================
# BƯỚC 3: Annotate (SnpEff/SnpSift)
# =========================
#echo "🧾 Annotate VCF: ${HAPLO_VCF}"

# Đảm bảo resources có index
#for v in "${GNOMAD_VCF}" "${CLINVAR_VCF}" "${THOUSANDG_VCF}"; do
#  [ -f "$v" ] || { echo "❌ Thiếu resource: $v"; exit 1; }
#  [[ "$v" =~ \.vcf\.gz$ ]] && { [ -f "${v}.tbi" ] || tabix -f -p vcf "$v"; }
#done

# Đảm bảo DB SnpEff
#if [ ! -d "${SNPEFF_HOME}/data/${SNPEFF_DB}" ]; then
  #echo "ℹ️ Tải snpEff DB: ${SNPEFF_DB}"
#  java ${JAVA_OPTS_SNPEFF} -jar "${SNPEFF_JAR}" download "${SNPEFF_DB}" -c "${SNPEFF_CONFIG}"
#fi

#echo "🔬 [1] snpEff → ${OUT_SNPEFF}"
#java ${JAVA_OPTS_SNPEFF} -jar "${SNPEFF_JAR}" ann -c "${SNPEFF_CONFIG}" -v "${SNPEFF_DB}" \
#  -stats "${SNPEFF_STATS}" \
#  "${HAPLO_VCF}" > "${OUT_SNPEFF}"
#
#echo "📊 [2] gnomAD (AF) → ${OUT_GNOMAD}"
#java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -id -info AF \
 # "${GNOMAD_VCF}" "${OUT_SNPEFF}" > "${OUT_GNOMAD}"

# Đổi AF -> GNOMAD_AF ngay sau gnomAD (nếu có bcftools)
#if command -v bcftools >/dev/null 2>&1; then
#  echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="Allele frequency from gnomAD">' > "${ANN_DIR}/gn_hdr.hdr"
#  bcftools annotate -h "${ANN_DIR}/gn_hdr.hdr" \
#    -c INFO/GNOMAD_AF:=INFO/AF -x INFO/AF \
 #   -O v -o "${OUT_GNOMAD}.tmp" "${OUT_GNOMAD}"
#  mv -f "${OUT_GNOMAD}.tmp" "${OUT_GNOMAD}"
 # rm -f "${ANN_DIR}/gn_hdr.hdr"
#else
#  echo "ℹ️ Không có bcftools → giữ nguyên INFO/AF từ gnomAD."
#fi

#echo "🧬 [3] ClinVar (CLNSIG, CLNDN) → ${OUT_CLINVAR}"
#java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate \
#  -info CLNSIG,CLNDN "${CLINVAR_VCF}" "${OUT_GNOMAD}" > "${OUT_CLINVAR}"

#echo "🌍 [4] 1000 Genomes (AF) → ${OUT_FINAL}"
#java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate \
#  -info AF "${THOUSANDG_VCF}" "${OUT_CLINVAR}" > "${OUT_FINAL}"

# Đổi AF -> KG_AF ngay sau 1000G (nếu có bcftools)
#if command -v bcftools >/dev/null 2>&1; then
#  echo '##INFO=<ID=KG_AF,Number=A,Type=Float,Description="Allele frequency from 1000 Genomes Project">' > "${ANN_DIR}/kg_hdr.hdr"
#  bcftools annotate -h "${ANN_DIR}/kg_hdr.hdr" \
#    -c INFO/KG_AF:=INFO/AF -x INFO/AF \
#    -O v -o "${OUT_FINAL}.tmp" "${OUT_FINAL}"
#  mv -f "${OUT_FINAL}.tmp" "${OUT_FINAL}"
#  rm -f "${ANN_DIR}/kg_hdr.hdr"
#else
#  echo "ℹ️ Không có bcftools → giữ nguyên INFO/AF cho 1000G."
#fi

# (tuỳ chọn) nén + index final
#if command -v bgzip >/dev/null 2>&1; then
#  bgzip -f "${OUT_FINAL}"
#  tabix -f -p vcf "${OUT_FINAL}.gz"
#  echo "📦 Ghi: ${OUT_FINAL}.gz"
#else
#  echo "📄 Ghi: ${OUT_FINAL}"
#fi

# =========================
# BƯỚC 4: DỌN RÁC Annotate
# =========================
#if [ "${CLEAN_INTERMEDIATE}" = "true" ]; then
#  echo "🧹 Dọn file trung gian annotate..."
#  rm -f "${OUT_SNPEFF}" "${OUT_GNOMAD}" "${OUT_CLINVAR}"
#  rm -f "${OUT_GNOMAD}.tmp" "${OUT_FINAL}.tmp"
#  rm -f "${HAPLO_GVCF}" "${HAPLO_GVCF}.tbi" 
#  rm -f "${ANN_DIR}/gn_hdr.hdr" "${ANN_DIR}/kg_hdr.hdr" 2>/dev/null || true
  # nếu đã nén final, dọn bản .vcf thường
#  if [ -f "${OUT_FINAL}.gz" ]; then
#    rm -f "${OUT_FINAL}"
#  fi
#  echo "✅ Đã dọn file trung gian annotate."
#else
#  echo "ℹ️ CLEAN_INTERMEDIATE=false → giữ lại file trung gian annotate."
#fi

#echo "✅ HOÀN TẤT cho mẫu ${SAMPLE_NAME}"
#echo "    • VCF (GATK): ${HAPLO_VCF}"
#echo "    • VCF đã annotate: ${OUT_FINAL}$( [ -f "${OUT_FINAL}.gz" ] && echo '.gz' )")

