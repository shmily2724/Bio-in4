
SAMPLE=$1

# Tạo biến cho đường dẫn input/output
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/${SAMPLE}"
RAW_DIR="${PROJECT_DIR}/raw_data"
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
mkdir -p "$SAMPLE_DIR"
mkdir -p "$TRIM_DIR"
# File input
READ1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
READ2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

# File output
TRIMMED_R1="${TRIM_DIR}/${SAMPLE}_1_paired.fastq.gz"
UNPAIRED_R1="${TRIM_DIR}/${SAMPLE}_1_unpaired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE}_2_paired.fastq.gz"
UNPAIRED_R2="${TRIM_DIR}/${SAMPLE}_2_unpaired.fastq.gz"
eval "$(conda shell.bash hook)"
conda activate BRCA
trimmomatic PE -threads 8 -phred33 
"${READ1}" "${READ2}" \
  "${TRIMMED_R1}" "${UNPAIRED_R1}" \
  "${TRIMMED_R2}" "${UNPAIRED_R2}" \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
