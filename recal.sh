SAMPLE_NAME=$1 
#input, output các biến đường dẫn
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
BWA_DIR="${SAMPLE_DIR}/Bwa alignments"
RECAL_DIR="${SAMPLE_DIR}/recal"
REF="/media/shmily/writable/BRCA_project/reference/reference.fa"
KNOWN_SNP="/media/shmily/writable/BRCA_project/reference/known_sites/Homo_sapiens_assembly38.dbsnp138.chr.vcf.gz"
KNOWN_INDEL="/media/shmily/writable/BRCA_project/reference/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"
1000_INDEL= "/media/shmily/writable/BRCA_project/reference/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
#input
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"
#ouput 
RECAL_DATA_TABLE="${RECAL_DIR}/${SAMPLE_NAME}_recal_data.table"
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"
#chạy conda
eval "$(conda shell.bash hook)" 
conda activate GATK 
# BaseRecalibrator 
gatk BaseRecalibrator \
 -I "${BAM_DEDUP}" \
 -R "${REF}" \
 -knownSites "${KNOWN_SNP}" \
 -knownSites "${KNOWN_INDEL}"\
 -knownSites "${1000_INDEL}" \
 -o "${RECAL_DATA_TABLE}"
 #ApplyBQSR 
 gatk ApplyBQSR \
   -R "${REF}" \
   -I "${BAM_DEDUP}" \
   --bqsr-recal-file "${RECAL_DATA_TABLE}" \
   -O "${RECAL_BAM}"
  samtools index "${RECAL_BAM}"
