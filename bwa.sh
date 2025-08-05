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
# Sort BAM với picard 
SORTED_BAM="${BWA_DIR}/${SAMPLE_NAME}_sorted.bam"
conda activate BRCA
PICARD_JAR="/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar" #đường dẫn chuẩn của picard
# Gọi Picard bằng lệnh java trực tiếp với -Djava.awt.headless=true
java -Djava.awt.headless=true -jar "${PICARD_JAR}" 
    SortSam \ 
    I="${BAM}"\
    O="${SORTED_BAM}" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=2000000
#Mark duplicate
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"
DEDUP_METRIX= "${BWA_DIR}/${SAMPLE_NAME}_dedup_metrics.txt"
java -Djava.awt.headless=true -jar "${PICARD_JAR}" 
  MarkDuplicates \
  MAX_RECORDS_IN_RAM=2000000 \
  VALIDATION_STRINGENCY=SILENT \
  M=$dedup_metrics \
  I="${SORTED_BAM}" \
  O="${BAM_DEDUP}" 
    
