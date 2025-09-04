#!/usr/bin/env bash
# ===================================================================
# == So sánh GATK vs DeepVariant bằng bcftools isec (đã sửa lỗi)  ==
# ==  - Đồng bộ 'chr' giữa GATK & DV                              ==
# ==  - Xuất only_gatk / only_dv / common + summary               ==
# ==  - Đếm đúng sau khi nén bgzip                                ==
# ===================================================================
set -euo pipefail

# ===== KÍCH HOẠT CONDA (đồng bộ với script annotate) =====
eval "$(conda shell.bash hook)"
conda activate BCF

echo "[DEBUG] PATH=$PATH"
echo "[DEBUG] bcftools exe: $(command -v bcftools)"
bcftools --version | head -n1

# ========================== CẤU HÌNH ==========================
PROJECT_DIR="${PROJECT_DIR:-/media/shmily/writable/BRCA_project}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_DIR}/isec_out}"   # nơi chứa output isec
THREADS="${THREADS:-4}"

# Cách 'collapse' của bcftools isec: snps|indels|both|none (mặc định: both)
COLLAPSE="${COLLAPSE:-both}"

# Bật/tắt tạo file thống kê tổng hợp
MAKE_SUMMARY="${MAKE_SUMMARY:-true}"

# ======================== HÀM TIỆN ÍCH ========================
TS() { date '+%Y-%m-%d %H:%M:%S'; }
need_cmd()  { command -v "$1" >/dev/null 2>&1 || { echo "❌ Thiếu tool: $1"; exit 1; }; }

validate_collapse() {
  case "$COLLAPSE" in
    snps|indels|both|none) : ;;
    *) echo "❌ COLLAPSE phải là: snps|indels|both|none (hiện: $COLLAPSE)"; exit 1;;
  esac
}

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

# Đổi tiền tố chr cho 1 file để khớp phong cách của 1 file template
# $1=sample_file  $2=template_file  $3=out_file(.vcf.gz)
harmonize_to_template() {
  local sample="$1" template="$2" out="$3"
  if has_chr "$sample" && ! has_chr "$template"; then
    echo "[$(TS)] 🔁 Bỏ 'chr' để khớp template ($(basename "$template"))"
    bcftools annotate --rename-chrs "$MKMAP_RMCHR" -Oz -o "$out" "$sample"
    tabix -f -p vcf "$out"
  elif ! has_chr "$sample" && has_chr "$template"; then
    echo "[$(TS)] 🔁 Thêm 'chr' để khớp template ($(basename "$template"))"
    bcftools annotate --rename-chrs "$MKMAP_ADDCHR" -Oz -o "$out" "$sample"
    tabix -f -p vcf "$out"
  else
    if [[ "$sample" != "$out" ]]; then
      cp -f "$sample" "$out"
      tabix -f -p vcf "$out" || true
    fi
  fi
}

# Đảm bảo .vcf.gz có index; nếu là .vcf thì bgzip + index
ensure_bgzip_tabix() {
  local in="$1"
  if [[ "$in" == *.vcf.gz ]]; then
    [[ -f "${in}.tbi" ]] || { echo "[$(TS)] ▶ Index: tabix -p vcf $in"; tabix -p vcf "$in"; }
    echo "$in"
  elif [[ "$in" == *.vcf ]]; then
    echo "[$(TS)] ▶ Nén bgzip: $in"
    bgzip -@ "${THREADS}" -f "$in"
    tabix -f -p vcf "${in}.gz"
    echo "${in}.gz"
  else
    echo "[$(TS)] ❌ Không nhận diện phần mở rộng: $in" >&2
    return 1
  fi
}

# Đếm số biến thể (tự động hỗ trợ .vcf và .vcf.gz)
count_variants() {
  local v="$1"
  if [[ -s "$v" ]]; then
    bcftools view -H "$v" | wc -l | awk '{print $1}'
  elif [[ -s "${v}.gz" ]]; then
    bcftools view -H "${v}.gz" | wc -l | awk '{print $1}'
  else
    echo 0
  fi
}

# Tìm file GATK & DV đã annotate theo quy ước script annotate
# Trả về GATK_VCF, DV_VCF (có thể rỗng nếu thiếu)
find_pair_for_sample() {
  local sample="$1" snpeff_dir="$2"
  GATK_VCF=""
  DV_VCF=""
  local g="${snpeff_dir}/${sample}_final_annotated.vcf"
  local d="${snpeff_dir}/${sample}_dv_final_annotated.vcf"
  if   [[ -f "${g}.gz" ]]; then GATK_VCF="${g}.gz"
  elif [[ -f "${g}"    ]]; then GATK_VCF="${g}"
  fi
  if   [[ -f "${d}.gz" ]]; then DV_VCF="${d}.gz"
  elif [[ -f "${d}"    ]]; then DV_VCF="${d}"
  fi
}

# Chạy isec cho 1 sample
run_isec_one_sample() {
  local sample="$1" gatk_in="$2" dv_in="$3" out_dir="$4"
  mkdir -p "$out_dir"

  # Bảo đảm nén + index
  gatk_in="$(ensure_bgzip_tabix "$gatk_in")"
  dv_in="$(ensure_bgzip_tabix "$dv_in")"

  # Đồng bộ kiểu 'chr' giữa 2 file (dùng GATK làm template)
  local gatk_harm="${out_dir}/${sample}.gatk.harm.vcf.gz"
  local dv_harm="${out_dir}/${sample}.dv.harm.vcf.gz"
  harmonize_to_template "$gatk_in" "$gatk_in" "$gatk_harm"   # copy + index nếu đã khớp
  harmonize_to_template "$dv_in"   "$gatk_in" "$dv_harm"

  echo "[$(TS)] ▶ bcftools isec: $sample"
  echo "                 • GATK: $gatk_harm"
  echo "                 • DeepVariant: $dv_harm"
  echo "                 • Collapse: $COLLAPSE"
  echo "                 • Out: $out_dir"

  # Thư mục thô để nhận 0000/0001/0002
  local raw_dir="${out_dir}/raw"
  mkdir -p "$raw_dir"

  # Thứ tự file quan trọng: 0000 -> chỉ file1 (GATK), 0001 -> chỉ file2 (DV), 0002 -> chung
  bcftools isec -c "$COLLAPSE" -p "$raw_dir" "$gatk_harm" "$dv_harm"

  # Đổi tên dễ hiểu
  local only_gatk="${out_dir}/${sample}.only_gatk.vcf"
  local only_dv="${out_dir}/${sample}.only_dv.vcf"
  local common="${out_dir}/${sample}.common.vcf"

  mv -f "${raw_dir}/0000.vcf" "$only_gatk"
  mv -f "${raw_dir}/0001.vcf" "$only_dv"
  mv -f "${raw_dir}/0002.vcf" "$common"

  # Nén + index cho output (cho tiện mở)
  bgzip -@ "${THREADS}" -f "$only_gatk"; tabix -f -p vcf "${only_gatk}.gz"
  bgzip -@ "${THREADS}" -f "$only_dv";   tabix -f -p vcf "${only_dv}.gz"
  bgzip -@ "${THREADS}" -f "$common";    tabix -f -p vcf "${common}.gz"

  # Tạo summary riêng (đếm **trên file .gz** để tránh 0)
  if [[ "$MAKE_SUMMARY" == "true" ]]; then
    local total_gatk total_dv c_og c_od c_cm
    total_gatk=$(count_variants "$gatk_harm")
    total_dv=$(count_variants "$dv_harm")
    c_og=$(count_variants "${only_gatk}.gz")
    c_od=$(count_variants "${only_dv}.gz")
    c_cm=$(count_variants "${common}.gz")

    echo -e "Sample\tTotal_GATK\tTotal_DV\tOnly_GATK\tOnly_DV\tCommon" > "${out_dir}/${sample}.summary.tsv"
    echo -e "${sample}\t${total_gatk}\t${total_dv}\t${c_og}\t${c_od}\t${c_cm}" >> "${out_dir}/${sample}.summary.tsv"

    # Ghi vào file tạm để tổng hợp cuối cùng
    echo -e "${sample}\t${total_gatk}\t${total_dv}\t${c_og}\t${c_od}\t${c_cm}" >> "${OUT_ROOT}/_project_summary.tmp"
  fi

  # Xoá raw tạm
  rm -rf "$raw_dir"
}

# =================== KIỂM TRA TIỀN ĐIỀU KIỆN ===================
for c in bcftools tabix bgzip; do need_cmd "$c"; done
validate_collapse
mkdir -p "$OUT_ROOT"

# ======================= MAIN (QUÉT/RUN) =======================
main() {
  # Chuẩn bị file tổng hợp
  if [[ "$MAKE_SUMMARY" == "true" ]]; then
    : > "${OUT_ROOT}/_project_summary.tmp"
  fi

  # Lấy danh sách sample
  local samples=()
  if (( $# > 0 )); then
    samples=("$@")
  else
    shopt -s nullglob
    for sdir in "${RESULTS_DIR}"/*; do
      [[ -d "$sdir" ]] || continue
      samples+=("$(basename "$sdir")")
    done
    shopt -u nullglob
  fi

  if (( ${#samples[@]} == 0 )); then
    echo "[$(TS)] ❌ Không tìm thấy sample nào trong: ${RESULTS_DIR}"
    exit 1
  fi

  local n_processed=0 n_skipped=0
  for S in "${samples[@]}"; do
    local S_DIR="${RESULTS_DIR}/${S}"
    local SNPEFF_DIR="${S_DIR}/snpeff"

    if [[ ! -d "$SNPEFF_DIR" ]]; then
      echo "[$(TS)] ⚠ Bỏ qua ${S}: không có thư mục snpeff/"
      ((n_skipped++)) || true
      continue
    fi

    find_pair_for_sample "$S" "$SNPEFF_DIR"
    if [[ -z "${GATK_VCF}" || -z "${DV_VCF}" ]]; then
      echo "[$(TS)] ⚠ Bỏ qua ${S}: thiếu ${S}_final_annotated(.vcf/.vcf.gz) hoặc ${S}_dv_final_annotated(.vcf/.vcf.gz)"
      ((n_skipped++)) || true
      continue
    fi

    local OUT_DIR="${OUT_ROOT}/${S}"
    run_isec_one_sample "$S" "$GATK_VCF" "$DV_VCF" "$OUT_DIR"
    ((n_processed++)) || true
  done

  if [[ "$MAKE_SUMMARY" == "true" ]]; then
    local summary="${OUT_ROOT}/project_isec_summary.tsv"
    echo -e "Sample\tTotal_GATK\tTotal_DV\tOnly_GATK\tOnly_DV\tCommon" > "$summary"
    cat "${OUT_ROOT}/_project_summary.tmp" >> "$summary"
    rm -f "${OUT_ROOT}/_project_summary.tmp"
    echo "[$(TS)] ✅ Tổng hợp: $summary"
  fi

  echo "[$(TS)] ✅ Hoàn tất. Processed: ${n_processed}, Skipped: ${n_skipped}"
}

main "$@"

# Dọn map tạm
rm -f "$MKMAP_ADDCHR" "$MKMAP_RMCHR"
