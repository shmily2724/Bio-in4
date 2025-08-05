SAMPLE_NAME=$1 
#input, output các biến đường dẫn
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
RECAL_DIR="${SAMPLE_DIR}/recal"
HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"
REF="/media/shmily/writable/BRCA_project/reference/reference.fa"
#input
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"
#ouput
HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"
mkdir -p ${HAPLO_DIR}
#chạy conda
eval "$(conda shell.bash hook)" 
conda activate GATK 
#GATK Haplotypecaller default 
gatk HaplotypeCaller \
 -R "${REF}" \
 -I "${RECAL_BAM}" \
 -O "${HAPLO_GVCF}" \
 -ERC GVCF
# GVCF -> VCF
conda activate GATK
gatk GenotypeGVCFs \
 -R "${REF}" \
 -V "${HAPLO_GVCF}" \
 -O "${HAPLO_VCF}"
# index GATK VCF
tabix -p vcf ${HAPLO_VCF}
