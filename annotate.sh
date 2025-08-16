#!/usr/bin/env bash
set -euo pipefail

# ========================== CẤU HÌNH ==========================
PROJECT_DIR="/media/shmily/writable/BRCA_project"

# Tham chiếu & resources
REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
GNOMAD_VCF="${PROJECT_DIR}/reference/resources/gnomad.v4.1.panel.merged.vcf.gz"
CLINVAR_VCF="${PROJECT_DIR}/reference/resources/clinvar_20250810.vcf.gz"
ONEKG_VCF="${PROJECT_DIR}/reference/resources/1000g.panel.merged.vcf.gz"

# SnpEff / SnpSift
SNPEFF_HOME="${PROJECT_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
SNPEFF_DB="GRCh38.86"

# Tuỳ chọn
CLEAN_INTERMEDIATE="${CLEAN_INTERMEDIATE:-true}"  # xóa file tạm sau khi chạy xong
GZIP_FINAL="${GZIP_FINAL:-true}"                  # nén file cuối
THREADS="${THREADS:-4}"
ATOMIZE="${ATOMIZE:-true}"                        # true|false – bật/tắt bước atomize

# ======================== HÀM TIỆN ÍCH ========================
TS() { date '+%Y-%m-%d %H:%M:%S'; }

need_cmd()  { command -v "$1" >/dev/null 2>&1 || { echo "❌ Thiếu tool: $1"; exit 1; }; }
need_file() { [[ -f "$1" ]] || { echo "❌ Thiếu file: $1"; exit 1; }; }

# header có 'ID=chr'?
has_chr() {
  bcftools view -h "$1" | grep -m1 '^##contig' | grep -q 'ID=chr' && return 0 || return 1
}

# Map rename chr <-> non-chr (tạo 1 lần)
MKMAP_ADDCHR="$(mktemp)"; cat > "$MKMAP_ADDCHR" <<'EOF'
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
X chrX
Y chrY
MT chrM
EOF
MKMAP_RMCHR="$(mktemp)"; awk '{print $2"\t"$1}' "$MKMAP_ADDCHR" > "$MKMAP_RMCHR"

# Đổi tiền tố chr của sample để khớp resource (giữa “chr1” và “1”)
# $1=sample.vcf.gz  $2=resource.vcf.gz  $3=out.vcf.gz
harmonize_to_resource() {
  local sample="$1" resource="$2" out="$3"
  if has_chr "$sample" && ! has_chr "$resource"; then
    echo "[$(TS)] 🔁 Bỏ 'chr' để khớp $(basename "$resource")"
    bcftools annotate --rename-chrs "$MKMAP_RMCHR" -O z -o "$out" "$sample"
    tabix -f "$out"
  elif ! has_chr "$sample" && has_chr "$resource"; then
    echo "[$(TS)] 🔁 Thêm 'chr' để khớp $(basename "$resource")"
    bcftools annotate --rename-chrs "$MKMAP_ADDCHR" -O z -o "$out" "$sample"
    tabix -f "$out"
  else
    if [[ "$sample" != "$out" ]]; then
      cp -f "$sample" "$out"
      tabix -f "$out" || true
    fi
  fi
}

# bcftools có hỗ trợ --atomize?
supports_bcftools_atomize() { bcftools norm -h 2>&1 | grep -q -- '--atomize'; }

# Atomize (tách MNP/complex thành primitives) sau khi đã norm
# $1=in.vcf.gz  $2=out.vcf.gz
atomize_after_norm() {
  local in_gz="$1" out_gz="$2"
  if [[ "$ATOMIZE" != "true" ]]; then
    echo "[$(TS)] ⏭️  Bỏ qua atomize (ATOMIZE=false)"; cp -f "$in_gz" "$out_gz"; tabix -f "$out_gz" || true; return
  fi
  if supports_bcftools_atomize; then
    echo "[$(TS)] 🧩 Atomize bằng bcftools --atomize..."
    bcftools norm --atomize -f "$REF" -O z -o "$out_gz" "$in_gz"
    tabix -f "$out_gz"
  elif command -v vt >/dev/null 2>&1; then
    echo "[$(TS)] 🧩 Atomize bằng vt decompose -s..."
    local tmp_vcf="${out_gz%.gz}"
    bcftools view -Ov -o "$tmp_vcf" "$in_gz"
    vt decompose -s "$tmp_vcf" -o "${tmp_vcf%.vcf}.atom.vcf"
    bgzip -f "${tmp_vcf%.vcf}.atom.vcf"; tabix -f -p vcf "${tmp_vcf%.vcf}.atom.vcf.gz"
    mv -f "${tmp_vcf%.vcf}.atom.vcf.gz" "$out_gz"
    mv -f "${tmp_vcf%.vcf}.atom.vcf.gz.tbi" "${out_gz}.tbi"
    rm -f "$tmp_vcf"
  else
    echo "[$(TS)] ⚠️ Không có bcftools --atomize hoặc vt → tiếp tục KHÔNG atomize."
    cp -f "$in_gz" "$out_gz"; tabix -f "$out_gz" || true
  fi
}

# Đổi INFO/AF thành tag mới (GNOMAD_AF hoặc KG_AF) rồi xoá AF gốc
# $1=input.vcf(.gz)  $2=new_tag  $3=out.vcf(.gz|.vcf)
rename_AF_to_new_tag() {
  local in="$1" newtag="$2" out="$3"
  local hdr tmp_tsv tsv_gz
  hdr="$(mktemp)"
  echo "##INFO=<ID=${newtag},Number=A,Type=Float,Description=\"Allele frequency from ${newtag}\">" > "$hdr"
  tmp_tsv="$(mktemp)"
  bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$in" > "$tmp_tsv"
  tsv_gz="${tmp_tsv}.gz"; bgzip -f -c "$tmp_tsv" > "$tsv_gz" && tabix -f -s 1 -b 2 -e 2 "$tsv_gz"
  bcftools annotate \
    -a "$tsv_gz" -c CHROM,POS,REF,ALT,INFO/"$newtag" \
    -h "$hdr" -x INFO/AF \
    -O v -o "$out" "$in"
  rm -f "$hdr" "$tmp_tsv" "$tsv_gz" "${tsv_gz}.tbi"
}

# Cần REF.fai cho bcftools norm
ensure_ref_index() {
  if [[ ! -f "${REF}.fai" ]]; then
    if command -v samtools >/dev/null 2>&1; then
      echo "[$(TS)] 🧩 Tạo index FASTA tham chiếu..."
      samtools faidx "$REF"
    else
      echo "❌ Thiếu ${REF}.fai và không có samtools để tạo. Cài samtools rồi chạy lại."
      exit 1
    fi
  fi
}

# =================== KIỂM TRA TIỀN ĐIỀU KIỆN ===================
for c in bcftools tabix bgzip java; do need_cmd "$c"; done
need_file "$REF"
need_file "$GNOMAD_VCF"
need_file "$CLINVAR_VCF"
need_file "$ONEKG_VCF"
need_file "$SNPEFF_JAR"
need_file "$SNPSIFT_JAR"

# Tải DB SnpEff nếu chưa có
if [[ ! -d "${SNPEFF_HOME}/data/${SNPEFF_DB}" ]]; then
  echo "[$(TS)] ⚠️ Chưa có DB SnpEff ${SNPEFF_DB} → tải về..."
  java -jar "$SNPEFF_JAR" download "$SNPEFF_DB" -c "$SNPEFF_CONFIG"
fi

ensure_ref_index

# ======================= THAM SỐ DÒNG LỆNH =======================
if [[ $# -lt 1 ]]; then
  cat <<EOF
Cách dùng:
  $(basename "$0") [--which=gatk|dv|both] SAMPLE1 [SAMPLE2 ...]
Mặc định --which=both
EOF
  exit 1
fi

WHICH="both"
if [[ "${1-}" =~ ^--which= ]]; then
  WHICH="${1#--which=}"
  shift
fi

# ============================ CORE ============================
# Annotate một VCF của 1 sample
# $1=sample_id  $2=tag("gatk"|"dv")  $3=input.vcf.gz
annotate_one() {
  local S="$1" TAG="$2" VCF_IN="$3"
  local SAMPLE_DIR="${PROJECT_DIR}/results/${S}"
  local ANN_DIR="${SAMPLE_DIR}/snpeff"
  mkdir -p "$ANN_DIR"

  local SUF=""; [[ "$TAG" == "dv" ]] && SUF="_dv"
  echo; echo "[$(TS)] ▶️ Xử lý: ${S} (${TAG})"

  # 1) Normalize + split multi-allelic
  local NORM="${ANN_DIR}/${S}${SUF}.norm.vcf.gz"
  echo "[$(TS)] 🔧 Normalize + split multi-allelic..."
  bcftools norm -m -both -f "$REF" -O z -o "$NORM" "$VCF_IN"
  tabix -f "$NORM"

  # 1b) Atomize primitives (tách MNP/complex)
  local ATOM="${ANN_DIR}/${S}${SUF}.atom.vcf.gz"
  atomize_after_norm "$NORM" "$ATOM"

  # 2) Đồng bộ 'chr' với gnomAD
  local HARM="${ANN_DIR}/${S}${SUF}.harm.vcf.gz"
  harmonize_to_resource "$ATOM" "$GNOMAD_VCF" "$HARM"

  # 3) SnpEff (tạo trường ANN)
  local SNPEFF_VCF="${ANN_DIR}/${S}${SUF}_snpeff.vcf"
  echo "[$(TS)] 🔬 SnpEff..."
  java -Xmx6g -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
       -stats "${ANN_DIR}/${S}${SUF}" \
       "$HARM" > "$SNPEFF_VCF"

  # 4) gnomAD annotate (AF) → GNOMAD_AF
  local GNOMAD_STEP="${ANN_DIR}/${S}${SUF}_gnomad.vcf"
  echo "[$(TS)] 📊 gnomAD annotate..."
  java -Xmx4g -jar "$SNPSIFT_JAR" annotate -info AF "$GNOMAD_VCF" "$SNPEFF_VCF" > "$GNOMAD_STEP"
  local GNOMAD_RENAMED="${ANN_DIR}/${S}${SUF}_gnomad_renamed.vcf"
  echo "[$(TS)] 📝 Đổi AF → GNOMAD_AF..."
  rename_AF_to_new_tag "$GNOMAD_STEP" "GNOMAD_AF" "$GNOMAD_RENAMED"

  # 5) ClinVar annotate (CLNSIG, CLNDN)
  local CLINVAR_STEP="${ANN_DIR}/${S}${SUF}_clinvar.vcf"
  echo "[$(TS)] 🧬 ClinVar annotate..."
  java -Xmx4g -jar "$SNPSIFT_JAR" annotate -info CLNSIG,CLNDN "$CLINVAR_VCF" "$GNOMAD_RENAMED" > "$CLINVAR_STEP"

  # 6) 1000G annotate (AF) → KG_AF
  local ONEKG_STEP="${ANN_DIR}/${S}${SUF}_1kg.vcf"
  echo "[$(TS)] 🌍 1000G annotate..."
  java -Xmx4g -jar "$SNPSIFT_JAR" annotate -info AF "$ONEKG_VCF" "$CLINVAR_STEP" > "$ONEKG_STEP"

  # 7) Đổi AF → KG_AF & hoàn tất
  local FINAL_VCF="${ANN_DIR}/${S}${SUF}_final_annotated.vcf"
  echo "[$(TS)] 📝 Đổi AF → KG_AF..."
  rename_AF_to_new_tag "$ONEKG_STEP" "KG_AF" "$FINAL_VCF"

  # 8) Nén + index (tuỳ chọn)
  if [[ "$GZIP_FINAL" == "true" ]]; then
    bgzip -f "$FINAL_VCF"
    tabix -f -p vcf "${FINAL_VCF}.gz"
    echo "[$(TS)] ✅ Final: ${FINAL_VCF}.gz"
  else
    echo "[$(TS)] ✅ Final: ${FINAL_VCF}"
  fi

  # 9) Dọn file tạm (tuỳ chọn)
  if [[ "$CLEAN_INTERMEDIATE" == "true" ]]; then
    echo "[$(TS)] 🧹 Dọn file tạm..."
    rm -f "$NORM" "${NORM}.tbi" "$ATOM" "${ATOM}.tbi" \
          "$HARM" "${HARM}.tbi" \
          "$SNPEFF_VCF" "$GNOMAD_STEP" "$GNOMAD_RENAMED" \
          "$CLINVAR_STEP" "$ONEKG_STEP"
  fi
}

# ======================= CHẠY THEO SAMPLE =======================
if [[ $# -lt 1 ]]; then echo "❌ Thiếu SAMPLE_ID"; exit 1; fi
for S in "$@"; do
  SAMPLE_DIR="${PROJECT_DIR}/results/${S}"
  VC_GATK="${SAMPLE_DIR}/haplotypecaller/${S}_gatk.vcf.gz"
  VC_DV="${SAMPLE_DIR}/Deepvariants/${S}_deepvariant.vcf.gz"
  [[ -d "$SAMPLE_DIR" ]] || { echo "❌ Không thấy thư mục sample: $SAMPLE_DIR"; continue; }

  case "$WHICH" in
    gatk)
      need_file "$VC_GATK"
      annotate_one "$S" "gatk" "$VC_GATK"
      ;;
    dv)
      need_file "$VC_DV"
      annotate_one "$S" "dv" "$VC_DV"
      ;;
    both)
      if [[ -f "$VC_GATK" ]]; then annotate_one "$S" "gatk" "$VC_GATK"; else echo "⚠️ Bỏ GATK: thiếu $VC_GATK"; fi
      if [[ -f "$VC_DV"   ]]; then annotate_one "$S" "dv"   "$VC_DV";   else echo "⚠️ Bỏ DV:   thiếu $VC_DV";   fi
      ;;
    *)
      echo "❌ --which phải là gatk|dv|both"; exit 1;;
  esac
done

# Dọn map tạm
rm -f "$MKMAP_ADDCHR" "$MKMAP_RMCHR"
