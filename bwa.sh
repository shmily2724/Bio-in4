SAMPLE_NAME=$1
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/${SAMPLE}"

# File input
READ1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
READ2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

# File output
TRIMMED_R1="${TRIM_DIR}/${SAMPLE}_1_paired.fastq.gz"
UNPAIRED_R1="${TRIM_DIR}/${SAMPLE}_1_unpaired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE}_2_paired.fastq.gz"
UNPAIRED_R2="${TRIM_DIR}/${SAMPLE}_2_unpaired.fastq.gz"
RG_ID="${SAMPLE_NAME}_RG"
PLATFORM="Illumina" # Hoặc thay thế bằng nền tảng giải trình tự thực tế của bạn (ví dụ: PacBio, OxfordNanopore)
LIBRARY_ID="Lib1" # Thay đổi nếu bạn có nhiều thư viện cho cùng một mẫu
PLATFORM_UNIT="${SAMPLE_NAME}_${PLATFORM}_${LIBRARY_ID}" # Tạo một ID đơn vị nền tảng duy nhất
eval "$(conda shell.bash hook)"
conda activate BRCA
bwa mem -Y \
  -K 100000000 \
  -t 16 \
  -R "$rg_string" \
  "$reference_fasta_file" \
  "$fastq_file1" \
  "$fastq_file2" | \
samtools view -Shb -o "$bam_file" -
