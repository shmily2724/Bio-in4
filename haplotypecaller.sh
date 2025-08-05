SAMPLE_NAME=$1 
#input, output các biến đường dẫn
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
RECAL_DIR="${SAMPLE_DIR}/recal"
HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"
ANN_DIR="${SAMPLE_DIR}/snpeff"
REF="/media/shmily/writable/BRCA_project/reference/reference.fa"
#input
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"
#ouput
HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"
mkdir -p ${HAPLO_DIR}
mkdir -p ${ANN_DIR}
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

#Annotation with SnpEff 
conda activate GATK 
SNPEFF_JAR="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.jar"
SNPEFF_CONFIG="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.config"
SNPEFF_DB="GRCh38.86"
#input
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"
#ouput
ANN_VCF="${ANN_DIR}/${SAMPLE_NAME}_gatk_annotated.vcf"
#chạy snpeff
java -Xmx6g -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
"${HAPLO_VCF}" >"${ANN_VCF}" \
