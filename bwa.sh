SAMPLE_NAME=$1
#tạo các biến đường dẫn
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
BWA_DIR="${SAMPLE_DIR}/Bwa alignments"
REF="/media/shmily/writable/BRCA_project/reference/reference.fa"
# File input
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"
# File output
BAM="${BWA_DIR}/${SAMPLE_NAME}_aligned.bam"
mkdir-p "$BWA_DIR"
# Định nghĩa thông tin Read Group
# Các biến ${SAMPLE_NAME}, ${THREADS}, ${REF}, ${TRIMMED_R1}, ${TRIMMED_R2}, ${OUTDIR} 
# cần được định nghĩa trước đó trong script của bạn.
RG_ID="${SAMPLE_NAME}_RG"
PLATFORM="Illumina" # Hoặc thay thế bằng nền tảng giải trình tự thực tế của bạn (ví dụ: PacBio, OxfordNanopore)
LIBRARY_ID="Lib1" # Thay đổi nếu bạn có nhiều thư viện cho cùng một mẫu
PLATFORM_UNIT="${SAMPLE_NAME}_${PLATFORM}_${LIBRARY_ID}" # Tạo một ID đơn vị nền tảng duy nhất
eval "$(conda shell.bash hook)"
conda activate BRCA
bwa mem -Y \
  -K 100000000 \
  -t 16 \
  -R "@RG\tID:${RG_ID}\tSM:${SAMPLE_NAME}\tPL:${PLATFORM}\tLB:${LIBRARY_ID}\tPU:${PLATFORM_UNIT}" \
  "$REF" \
  "${TRIMMED_R1}" \
  "${TRIMMED_R2}" | \
samtools view -Shb -o "${BAM}" -
