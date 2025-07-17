#!/bin/bash

#=== Step 0: Logging utilities ===#
function log_step {
  STEP_NAME="$1"
  echo -e "\n\033[1;34m>>> ${STEP_NAME}...\033[0m"
}

function log_done {
  echo -e "\033[1;32m[Done]\033[0m"
}

function log_error {
  ERROR_MESSAGE="$1"
  echo -e "\033[1;31mERROR: ${ERROR_MESSAGE}\033[0m"
  exit 1
}

# Ngắt script ngay lập tức nếu bất kỳ lệnh nào thất bại
set -e

#=== Step 1: Check tools in Conda environments ===#
log_step "Step 1: Checking required tools in Conda environments"

function check_tool {
  TOOL_NAME="$1"
  ENV_NAME="$2"

  if ! conda run -n "$ENV_NAME" which "$TOOL_NAME" &>/dev/null; then
    echo -e "\033[1;31mERROR: $TOOL_NAME NOT FOUND in environment $ENV_NAME. Please install it.\033[0m"
    exit 1
  else
    echo -e "\033[1;32m✔ $TOOL_NAME found in environment $ENV_NAME\033[0m"
  fi
}

# Tools in BRCA env
check_tool trimmomatic BRCA
check_tool bwa BRCA
check_tool samtools BRCA
check_tool picard BRCA
check_tool mosdepth BRCA
check_tool fastqc BRCA
check_tool multiqc BRCA

# Tools in GATK env
check_tool gatk GATK
check_tool snpEff GATK

log_done

#=== Step 2: Input Configuration ===#
OUTDIR_NAME="$1" # Tên thư mục đầu ra, vd: ERR152
READ1_FILENAME="$2" # Tên file đọc 1, vd: ERR2228152_1.fastq.gz
READ2_FILENAME="$3" # Tên file đọc 2, vd: ERR2228152_2.fastq.gz

if [[ -z "$OUTDIR_NAME" || -z "$READ1_FILENAME" || -z "$READ2_FILENAME" ]]; then
  log_error "Missing arguments. Usage: bash pipeline_default.sh <OUTPUT_FOLDER_NAME> <READ1_FILENAME> <READ2_FILENAME>"
fi

# Trích xuất SAMPLE_NAME từ tên file đọc 1 (bỏ phần _1.fastq.gz)
# Giả định tên file đọc 1 luôn kết thúc bằng _1.fastq.gz
SAMPLE_NAME=$(basename "${READ1_FILENAME}" | sed 's/_1\.fastq\.gz$//')

READ1="raw_data/${READ1_FILENAME}"
READ2="raw_data/${READ2_FILENAME}"
REF="/media/shmily/writable/BRCA_project/reference/reference.fa"
KNOWN_SITES_DIR="/media/shmily/writable/BRCA_project/reference/known_sites"
OUTDIR="results/${OUTDIR_NAME}"

# Định nghĩa DeepVariant Docker Image Version
DEEPVARIANT_VERSION="1.9.0" # Hoặc phiên bản DeepVariant bạn đã pull

# Kiểm tra sự tồn tại của các file đầu vào
if [[ ! -f "$READ1" ]]; then
  log_error "Input Read1 file not found: ${READ1}. Please ensure it's in the 'raw_data' directory."
fi
if [[ ! -f "$READ2" ]]; then
  log_error "Input Read2 file not found: ${READ2}. Please ensure it's in the 'raw_data' directory."
fi
if [[ ! -f "$REF" ]]; then
  log_error "Reference file not found: ${REF}."
fi
if [[ ! -f "${KNOWN_SITES_DIR}/Mills_and_1000G_gold_standard.indels.hg38.chr_FIXED.vcf.gz" ]]; then
  log_error "Known sites VCF (Mills_and_1000G) not found at ${KNOWN_SITES_DIR}."
fi
if [[ ! -f "${KNOWN_SITES_DIR}/Homo_sapiens_assembly38.dbsnp138.chr_FIXED.vcf.gz" ]]; then
  log_error "Known sites VCF (dbsnp) not found at ${KNOWN_SITES_DIR}."
fi

mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/fastqc_raw
mkdir -p ${OUTDIR}/fastqc_trimmed
mkdir -p ${OUTDIR}/deepvariant_logs

#=== Bước tiền xử lý Conda hook (chỉ chạy 1 lần) ===#
eval "$(conda shell.bash hook)"

#=== Step 3: Quality Control on Raw Reads (FastQC) ===#
log_step "Step 3: Running FastQC on raw reads"
conda activate BRCA
fastqc ${READ1} -o ${OUTDIR}/fastqc_raw 2>&1 | tee -a ${OUTDIR}/fastqc_raw.log
fastqc ${READ2} -o ${OUTDIR}/fastqc_raw 2>&1 | tee -a ${OUTDIR}/fastqc_raw.log
log_done

#=== Step 4: Quality trimming ===#
log_step "Step 4: Quality trimming with Trimmomatic"
conda activate BRCA
# Đầu ra của Trimmomatic sẽ là _1_paired.fq.gz và _2_paired.fq.gz
TRIMMED_R1="${OUTDIR}/${SAMPLE_NAME}_1_paired.fq.gz"
TRIMMED_R2="${OUTDIR}/${SAMPLE_NAME}_2_paired.fq.gz"
UNPAIRED_R1="${OUTDIR}/${SAMPLE_NAME}_1_unpaired.fq.gz"
UNPAIRED_R2="${OUTDIR}/${SAMPLE_NAME}_2_unpaired.fq.gz"

trimmomatic PE -threads 8 \
  ${READ1} ${READ2} \
  ${TRIMMED_R1} ${UNPAIRED_R1} \
  ${TRIMMED_R2} ${UNPAIRED_R2} \
  SLIDINGWINDOW:4:20 MINLEN:36 \
  2>&1 | tee -a ${OUTDIR}/trimmomatic.log
log_done

#=== Step 5: Quality Control on Trimmed Reads (FastQC) ===#
log_step "Step 5: Running FastQC on trimmed reads"
fastqc ${TRIMMED_R1} -o ${OUTDIR}/fastqc_trimmed 2>&1 | tee -a ${OUTDIR}/fastqc_trimmed.log
fastqc ${TRIMMED_R2} -o ${OUTDIR}/fastqc_trimmed 2>&1 | tee -a ${OUTDIR}/fastqc_trimmed.log
log_done

#=== Step 6: Alignment ===#
log_step "Step 6: BWA MEM alignment"

# Định nghĩa thông tin Read Group
# Các biến ${SAMPLE_NAME}, ${THREADS}, ${REF}, ${TRIMMED_R1}, ${TRIMMED_R2}, ${OUTDIR} 
# cần được định nghĩa trước đó trong script của bạn.
RG_ID="${SAMPLE_NAME}_RG"
PLATFORM="Illumina" # Hoặc thay thế bằng nền tảng giải trình tự thực tế của bạn (ví dụ: PacBio, OxfordNanopore)
LIBRARY_ID="Lib1" # Thay đổi nếu bạn có nhiều thư viện cho cùng một mẫu
PLATFORM_UNIT="${SAMPLE_NAME}_${PLATFORM}_${LIBRARY_ID}" # Tạo một ID đơn vị nền tảng duy nhất

# Sửa lỗi 'rm: No such file or directory' bằng cách kiểm tra trước khi xóa
# Lệnh 'rm' sẽ không báo lỗi nếu file không tồn tại nhờ '|| true'
rm "${UNPAIRED_R1}" || true
rm "${UNPAIRED_R2}" || true

# Chạy bwa mem với tùy chọn -R để thêm thông tin Read Group vào header của file BAM.
# Chuyển hướng trực tiếp output sang BAM bằng cách pipe qua samtools view để tối ưu hiệu suất và không gian đĩa.
bwa mem -t "${THREADS}" -M \
    -R "@RG\tID:${RG_ID}\tSM:${SAMPLE_NAME}\tPL:${PLATFORM}\tLB:${LIBRARY_ID}\tPU:${PLATFORM_UNIT}" \
    "${REF}" "${TRIMMED_R1}" "${TRIMMED_R2}" | \
    samtools view -Sb - \
    > "${OUTDIR}/${SAMPLE_NAME}_aligned.bam" \
    2>&1 | tee -a "${OUTDIR}/bwa_mem.log"

log_done

#=== Step 7: Convert SAM to BAM + Sort ===#
log_step "Step 7: Sort BAM and index..."

# Định nghĩa các đường dẫn file
INPUT_ALIGNED_BAM="${OUTDIR}/${SAMPLE_NAME}_aligned.bam" # Đây là output từ Bước 6 đã sửa
OUTPUT_SORTED_BAM="${OUTDIR}/${SAMPLE_NAME}_sorted.bam"
OUTPUT_SORTED_BAI="${OUTPUT_SORTED_BAM}.bai" # Tên file chỉ mục

# Kiểm tra xem file BAM đầu vào từ Bước 6 có tồn tại và không rỗng không
if [ ! -s "${INPUT_ALIGNED_BAM}" ]; then
    echo "Lỗi: File BAM đã căn chỉnh (${INPUT_ALIGNED_BAM}) từ Bước 6 bị thiếu hoặc trống."
    log_done_error
    exit 1
fi

# Sắp xếp file BAM
echo "Sắp xếp file BAM: ${INPUT_ALIGNED_BAM} -> ${OUTPUT_SORTED_BAM}"
samtools sort -@ "${THREADS}" "${INPUT_ALIGNED_BAM}" -o "${OUTPUT_SORTED_BAM}" \
    2>&1 | tee -a "${OUTDIR}/samtools_sort.log"

# Kiểm tra xem lệnh sort có thành công không
if [ $? -ne 0 ]; then
    echo "Lỗi: samtools sort gặp sự cố. Vui lòng kiểm tra ${OUTDIR}/samtools_sort.log để biết chi tiết."
    log_done_error
    exit 1
fi

# Lập chỉ mục file BAM đã sắp xếp (quan trọng cho các công cụ như GATK)
echo "Lập chỉ mục file BAM đã sắp xếp: ${OUTPUT_SORTED_BAM}"
samtools index -@ "${THREADS}" "${OUTPUT_SORTED_BAM}" \
    2>&1 | tee -a "${OUTDIR}/samtools_index.log"

# Kiểm tra xem lệnh index có thành công không
if [ $? -ne 0 ]; then
    echo "Lỗi: samtools index gặp sự cố. Vui lòng kiểm tra ${OUTDIR}/samtools_index.log để biết chi tiết."
    log_done_error
    exit 1
fi

# Xóa file BAM không sắp xếp để tiết kiệm không gian
echo "Xóa file BAM không sắp xếp: ${INPUT_ALIGNED_BAM}"
rm "${INPUT_ALIGNED_BAM}" || true

log_done

#=== Step 8: Mark Duplicates ===#
log_step "Step 8: Mark duplicates with Picard"
conda activate BRCA

# GÁN TRỰC TIẾP ĐƯỜNG DẪN CHÍNH XÁC CỦA PICARD.JAR TỪ KẾT QUẢ FIND CỦA BẠN
PICARD_JAR="/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar" # <<< ĐÂY LÀ ĐƯỜNG DẪN CHÍNH XÁC BẠN TÌM THẤY

if [[ ! -f "$PICARD_JAR" ]]; then # Kiểm tra an toàn: đảm bảo file JAR tồn tại
    log_error "Picard JAR file not found at ${PICARD_JAR}. Please verify the path or reinstall Picard."
fi

# Gọi Picard bằng lệnh java trực tiếp với -Djava.awt.headless=true
java -Djava.awt.headless=true -jar "${PICARD_JAR}" MarkDuplicates \
  I=${OUTDIR}/${SAMPLE_NAME}_sorted.bam \
  O=${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  M=${OUTDIR}/${SAMPLE_NAME}_metrics.txt \
  2>&1 | tee -a ${OUTDIR}/picard_markduplicates.log
samtools index ${OUTDIR}/${SAMPLE_NAME}_marked.bam 2>&1 | tee -a ${OUTDIR}/samtools_index_marked.log
rm ${OUTDIR}/${SAMPLE_NAME}_sorted.bam || true
log_done

#=== Step 9: Base Recalibration ===#
log_step "Step 9: Base Quality Score Recalibration"
conda activate GATK
gatk BaseRecalibrator \
  -I ${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  -R ${REF} \
  --known-sites ${KNOWN_SITES_DIR}/Mills_and_1000G_gold_standard.indels.hg38.chr_FIXED.vcf.gz \
  --known-sites ${KNOWN_SITES_DIR}/Homo_sapiens_assembly38.dbsnp138.chr_FIXED.vcf.gz \
  -O ${OUTDIR}/recal_data.table \
  2>&1 | tee -a ${OUTDIR}/gatk_baserecalibrator.log

gatk ApplyBQSR \
  -R ${REF} \
  -I ${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  --bqsr-recal-file ${OUTDIR}/recal_data.table \
  -O ${OUTDIR}/${SAMPLE_NAME}_recal.bam \
  2>&1 | tee -a ${OUTDIR}/gatk_applybqsr.log

# Xóa các file trung gian sau BQSR
rm ${OUTDIR}/recal_data.table || true
rm ${OUTDIR}/${SAMPLE_NAME}_marked.bam || true
rm ${OUTDIR}/${SAMPLE_NAME}_marked.bam.bai || true
log_done

#=== Step 10: Variant Calling with DeepVariant ===#
log_step "Step 10: Variant Calling with DeepVariant (using Docker)"

# Kiểm tra xem Docker có đang chạy không.
# Sử dụng 'docker info' để kiểm tra trạng thái Docker Daemon.
if ! docker info &>/dev/null; then
    log_error "Docker is not running or accessible. Please start Docker."
    log_done_error
    exit 1
fi

# Đảm bảo các đường dẫn đầu vào/đầu ra là đường dẫn tuyệt đối cho Docker volume mounting.
# 'realpath' sẽ chuyển đổi đường dẫn tương đối thành tuyệt đối.
OUTDIR_ABS=$(realpath "${OUTDIR}")
REF_DIR_ABS=$(realpath "$(dirname "${REF}")") # Lấy đường dẫn tuyệt đối của thư mục chứa REF

# Tên file reference bên trong container (phải khớp với cách nó được gắn)
# dirname "${REF}" sẽ là "/reference_data" bên trong container.
# basename "${REF}" sẽ là tên file reference (ví dụ: reference.fa).
REF_IN_CONTAINER="/reference_data/$(basename "${REF}")"

# Chạy DeepVariant bằng Docker
docker run \
    -v "${OUTDIR_ABS}":"/output_data" \
    -v "${REF_DIR_ABS}":"/reference_data" \
    google/deepvariant:"${DEEPVARIANT_VERSION}" \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="${REF_IN_CONTAINER}" \
    --reads="/output_data/${SAMPLE_NAME}_recal.bam" \
    --output_vcf="/output_data/${SAMPLE_NAME}_deepvariant.vcf.gz" \
    --output_gvcf="/output_data/${SAMPLE_NAME}_deepvariant.g.vcf.gz" \
    --num_shards=$(nproc) \
    --logging_dir="/output_data/deepvariant_logs" \
    2>&1 | tee -a "${OUTDIR}/deepvariant.log"

# Kiểm tra mã thoát của lệnh Docker
if [ $? -ne 0 ]; then
    echo "Lỗi: DeepVariant gặp sự cố. Vui lòng kiểm tra ${OUTDIR}/deepvariant.log để biết chi tiết."
    log_done_error
    exit 1
fi

log_done

#=== Step 11: Index DeepVariant VCF ===#
log_step "Step 11: Indexing DeepVariant VCF"
tabix -p vcf ${OUTDIR}/${SAMPLE_NAME}_deepvariant.vcf.gz 2>&1 | tee -a ${OUTDIR}/tabix_deepvariant_index.log
log_done

# !!! QUAN TRỌNG: LỰA CHỌN VỀ XÓA FILE _recal.bam !!!
# Nếu bạn muốn GATK HaplotypeCaller (Bước 12) và Mosdepth (Bước 16) chạy trên file _recal.bam,
# bạn phải BỎ COMMENT (hoặc xóa) các dòng 'rm' dưới đây để giữ lại file.
# Nếu không, các bước đó sẽ bị lỗi vì thiếu file.
# rm ${OUTDIR}/${SAMPLE_NAME}_recal.bam || true
# rm ${OUTDIR}/${SAMPLE_NAME}_recal.bam.bai || true

#=== Step 12: Variant Calling with GATK HaplotypeCaller (tùy chọn) ===#
log_step "Step 12: Variant Calling with GATK HaplotypeCaller (Optional - Can be removed if DeepVariant is preferred)"
conda activate GATK
# Đảm bảo file _recal.bam tồn tại nếu bạn chạy bước này
gatk HaplotypeCaller \
  -R ${REF} \
  -I ${OUTDIR}/${SAMPLE_NAME}_recal.bam \
  -O ${OUTDIR}/${SAMPLE_NAME}_gatk.g.vcf.gz \
  -ERC GVCF \
  2>&1 | tee -a ${OUTDIR}/gatk_haplotypecaller.log
log_done  

#=== Step 13: Genotype GVCF (Nếu dùng GATK HaplotypeCaller) ===#
log_step "Step 13: Genotyping GATK GVCF to final VCF (Optional - If GATK HaplotypeCaller is used)"
conda activate GATK
gatk GenotypeGVCFs \
  -R ${REF} \
  -V ${OUTDIR}/${SAMPLE_NAME}_gatk.g.vcf.gz \
  -O ${OUTDIR}/${SAMPLE_NAME}_gatk.vcf.gz \
  2>&1 | tee -a ${OUTDIR}/gatk_genotypegvcfs.log
log_done

#=== Step 14: Index final GATK VCF (Nếu dùng GATK) ===#
log_step "Step 14: Indexing final GATK VCF (Optional - If GATK is used)"
tabix -p vcf ${OUTDIR}/${SAMPLE_NAME}_gatk.vcf.gz 2>&1 | tee -a ${OUTDIR}/tabix_gatk_index.log
log_done

#=== Step 15: Annotation with SnpEff ===#
log_step "Step 15: Annotating variants with SnpEff"
conda activate GATK
SNPEFF_CONFIG="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.config"
SNPEFF_DB="GRCh38.86"
# *** PHẦN ĐÃ THAY ĐỔI / BỔ SUNG: Cấp phát bộ nhớ cho SnpEff (-Xmx6g) và gọi trực tiếp Java ***

# Với 8GB RAM, 6GB là một lựa chọn hợp lý. Nếu gặp lỗi OutOfMemoryError nữa, hãy thử giảm xuống 5g

java -Xmx6g -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \

  "${OUTDIR}/${SAMPLE_NAME}_gatk.vcf.gz" > "${OUTDIR}/${SAMPLE_NAME}_gatk_annotated.vcf" \

  2>> "${OUTDIR}/snpeff_gatk_annotation.log" # Chuyển stderr vào file log



# Kiểm tra mã thoát của lệnh Java vừa rồi

if [ $? -ne 0 ]; then

    log_error "SnpEff annotation encountered an issue. Please check ${OUTDIR}/snpeff_gatk_annotation.log for details."

fi

log_done

#=== Step 16: Depth Analysis ===#
log_step "Step 16: Depth analysis with mosdepth"
conda activate BRCA
# Đảm bảo file _recal.bam tồn tại nếu bạn chạy bước này
mosdepth ${OUTDIR}/${SAMPLE_NAME} ${OUTDIR}/${SAMPLE_NAME}_recal.bam 2>&1 | tee -a ${OUTDIR}/mosdepth.log
log_done

#=== Step 17: MultiQC Report ===#
log_step "Step 17: Generating MultiQC report"
multiqc ${OUTDIR} -o ${OUTDIR} --filename multiqc_report.html 2>&1 | tee -a ${OUTDIR}/multiqc.log
log_done

#=== Finished ===#
echo -e "\n\033[1;35m>>> Pipeline completed for ${OUTDIR_NAME} <<<\033[0m"
echo -e "Check the MultiQC report: ${OUTDIR}/multiqc_report.html"
echo -e "DeepVariant VCF: ${OUTDIR}/${SAMPLE_NAME}_deepvariant.vcf.gz"
echo -e "DeepVariant Annotated VCF: ${OUTDIR}/${SAMPLE_NAME}_deepvariant_annotated.vcf"

