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
THREADS=8 # Đặt số luồng vào một biến để dễ thay đổi

# --- TẠO CÁC BIẾN ĐƯỜNG DẪN ---
PROJECT_DIR="/media/shmily/writable/BRCA_project"
ADAPTER_FILE="/home/shmily/miniconda/envs/BRCA/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

# Đã cập nhật biến SAMPLE_NAME
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
RAW_DIR="${PROJECT_DIR}/raw_data"
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"

# Thư mục cho kết quả FastQC
FASTQC_RAW_DIR="${SAMPLE_DIR}/fastqc_raw"
FASTQC_TRIM_DIR="${SAMPLE_DIR}/fastqc_trimmed"

# Tạo các thư mục output nếu chưa tồn tại
mkdir -p "${TRIM_DIR}"
mkdir -p "${FASTQC_RAW_DIR}"
mkdir -p "${FASTQC_TRIM_DIR}"

# --- Định nghĩa file Input và Output ---
# File input
READ1="${RAW_DIR}/${SAMPLE_NAME}_1.fastq.gz"
READ2="${RAW_DIR}/${SAMPLE_NAME}_2.fastq.gz"

# File output của Trimmomatic
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
UNPAIRED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_unpaired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"
UNPAIRED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_unpaired.fastq.gz"


# --- KÍCH HOẠT MÔI TRƯỜG CONDA ---
eval "$(conda shell.bash hook)"
conda activate BRCA


# --- BƯỚC 1: CHẠY FASTQC TRÊN DỮ LIỆU THÔ ---
echo "Bắt đầu chạy FastQC trên dữ liệu thô cho mẫu: ${SAMPLE_NAME}"

fastqc --threads ${THREADS} \
    -o "${FASTQC_RAW_DIR}" \
    "${READ1}" \
    "${READ2}"

echo "FastQC trên dữ liệu thô hoàn tất."


# --- BƯỚC 2: CHẠY TRIMMOMATIC ---
echo "Bắt đầu chạy Trimmomatic cho mẫu: ${SAMPLE_NAME}"

trimmomatic PE -threads ${THREADS} -phred33 \
    "${READ1}" "${READ2}" \
    "${TRIMMED_R1}" "${UNPAIRED_R1}" \
    "${TRIMMED_R2}" "${UNPAIRED_R2}" \
    ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

echo "Trimmomatic hoàn tất."


# --- BƯỚC 3: CHẠY FASTQC TRÊN DỮ LIỆU ĐÃ TRIM ---
echo "Bắt đầu chạy FastQC trên dữ liệu đã trim cho mẫu: ${SAMPLE_NAME}"

fastqc --threads ${THREADS} \
    -o "${FASTQC_TRIM_DIR}" \
    "${TRIMMED_R1}" \
    "${TRIMMED_R2}"

echo "FastQC trên dữ liệu đã trim hoàn tất."
echo "Pipeline trimming và QC hoàn tất cho mẫu: ${SAMPLE_NAME}!"
