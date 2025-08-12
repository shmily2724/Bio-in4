#!/bin/bash

# Tự động dừng script khi có lỗi
set -e

# Kiểm tra có ít nhất 1 tham số không
if [ "$#" -lt 1 ]; then
    echo "❌ Lỗi: Vui lòng cung cấp ít nhất một tên mẫu (SAMPLE_NAME)."
    echo "📌 Cách dùng: $0 <sample_name1> [sample_name2] [...]"
    exit 1
fi

# --- THIẾT LẬP CÁC BIẾN CỐ ĐỊNH ---
PROJECT_DIR="/media/shmily/writable/BRCA_project"
REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"
SNPEFF_HOME="${PROJECT_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
SNPEFF_DB="GRCh38.86"
SNPEFF_DATA_DIR="${SNPEFF_HOME}/data"
CLEAN_INTERMEDIATE=true  # ✅ Tuỳ chọn xóa file trung gian

# Kiểm tra database SnpEff đã tồn tại chưa
if [ ! -d "${SNPEFF_DATA_DIR}/${SNPEFF_DB}" ]; then
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] ⚠️ SnpEff DB ${SNPEFF_DB} chưa tồn tại. Đang tải về..."
    java -jar "$SNPEFF_JAR" download "$SNPEFF_DB" -c "$SNPEFF_CONFIG"
fi

# Lặp qua từng mẫu được truyền vào
for SAMPLE_NAME in "$@"; do
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] ▶️ Bắt đầu xử lý mẫu: ${SAMPLE_NAME}"

    SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
    HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"
    ANN_DIR="${SAMPLE_DIR}/snpeff"
    mkdir -p "${ANN_DIR}"

    # File input và output
    HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"
    ANN_VCF="${ANN_DIR}/${SAMPLE_NAME}_gatk_annotated.vcf"
    GNOMAD_ANN_VCF="${ANN_DIR}/${SAMPLE_NAME}_gnomad_annotated.vcf"
    CLINVAR_ANN_VCF="${ANN_DIR}/${SAMPLE_NAME}_clinvar_annotated.vcf"
    FINAL_ANN_VCF="${ANN_DIR}/${SAMPLE_NAME}_final_annotated.vcf"
    SNPEFF_STATS="${ANN_DIR}/${SAMPLE_NAME}"

    # Kiểm tra file đầu vào
    if [ ! -f "${HAPLO_VCF}" ]; then
        echo "❌ Lỗi: File VCF đầu vào không tồn tại: ${HAPLO_VCF}"
        continue
    fi

    # ✅ Chú giải với SnpEff
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] 🔬 Annotating với SnpEff..."
    java -Xmx6g -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
        -stats "${SNPEFF_STATS}" \
        "${HAPLO_VCF}" > "${ANN_VCF}"

    # ✅ Annotate với gnomAD
    GNOMAD_VCF="${PROJECT_DIR}/reference/resources/gnomad.v4.1.panel.merged.vcf.gz"
    if [ ! -f "${GNOMAD_VCF}" ]; then
        echo "❌ Lỗi: Không tìm thấy file gnomAD: ${GNOMAD_VCF}"
        continue
    fi
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] 📊 Annotating với gnomAD..."
    java -Xmx4g -jar "$SNPSIFT_JAR" annotate \
        -id -info AF "${GNOMAD_VCF}" "${ANN_VCF}" > "${GNOMAD_ANN_VCF}"

    # ✅ Annotate với ClinVar
    CLINVAR_VCF="${PROJECT_DIR}/reference/resources/clinvar_20250810.vcf.gz"
    if [ ! -f "${CLINVAR_VCF}" ]; then
        echo "❌ Lỗi: Không tìm thấy file ClinVar: ${CLINVAR_VCF}"
        continue
    fi
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] 🧬 Annotating với ClinVar..."
    java -Xmx4g -jar "$SNPSIFT_JAR" annotate \
        -info CLNSIG,CLNDN "${CLINVAR_VCF}" "${GNOMAD_ANN_VCF}" > "${CLINVAR_ANN_VCF}"

    # ✅ Annotate với 1000 Genomes
    THOUSANDG_VCF="${PROJECT_DIR}/reference/resources/1000g.panel.merged.vcf.gz"
    if [ ! -f "${THOUSANDG_VCF}" ]; then
        echo "❌ Lỗi: Không tìm thấy file 1000 Genomes: ${THOUSANDG_VCF}"
        continue
    fi
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] 🌍 Annotating với 1000 Genomes..."
    java -Xmx4g -jar "$SNPSIFT_JAR" annotate \
        -info AF "${THOUSANDG_VCF}" "${CLINVAR_ANN_VCF}" > "${FINAL_ANN_VCF}"

    # ✅ Xoá file trung gian nếu cần
    if [ "$CLEAN_INTERMEDIATE" = true ]; then
        echo "[`date '+%Y-%m-%d %H:%M:%S'`] 🧹 Xóa các file trung gian..."
        rm -f "$ANN_VCF" "$GNOMAD_ANN_VCF" "$CLINVAR_ANN_VCF"
    fi

    echo "[`date '+%Y-%m-%d %H:%M:%S'`] ✅ Hoàn tất pipeline cho mẫu: ${SAMPLE_NAME}"
done


