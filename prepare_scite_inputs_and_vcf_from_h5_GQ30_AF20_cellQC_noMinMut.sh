#!/bin/bash
set -uo pipefail

BASE="/Users/marumoatsushi/BCOR/SC"
H5DIR="${BASE}/resultstapestri/YoshidaDatah5"
SCITEDIR="${BASE}/SCITE"
SCRIPTDIR="${SCITEDIR}/scripts"

# New output directory to avoid overwriting the previous run
OUTBASE="${SCITEDIR}/output_DP10_GQ30_AF20_cellQC_noMinMut"

SCITE_OUT="${OUTBASE}"
VCF_OUT="${OUTBASE}/vcf"
LOG_OUT="${OUTBASE}/logs"

mkdir -p "${SCITE_OUT}" "${VCF_OUT}" "${LOG_OUT}"

# =========================
# QC settings based on the attached method text
# =========================
MIN_CELL_DEPTH=10
MIN_GENOTYPE_QUALITY=30
MIN_MUT_AF=20

# Variant-level filter: variants genotyped in fewer than 50% of cells are removed.
# In other words, keep variants with missing_fraction <= 0.50 after cell-level QC.
MAX_MISSING_FRAC=0.50

# Cell-level filter: cells with less than 50% of genotypes present are removed.
# In other words, keep cells with missing_fraction <= 0.50 after DP/GQ/AF rule.
MAX_CELL_MISSING_FRAC=0.50

# User requested: do NOT apply MIN_MUTATED_CELLS=5 at the h5-to-matrix step.
# This is set to 0 so it is not used. Mutation-frequency filtering can be done later on the final VCF.
MIN_MUTATED_CELLS=0

# Optional final VCF threshold from the method text: variants mutated in fewer than 1.0% of cells.
# This is NOT applied to the matrix by default. Use the optional command shown near the end if needed.
MIN_MUTATED_PERCENT_FOR_FINAL_VCF=1.0

# =========================
# Target h5 files
# =========================
H5_FILES=(
  "${H5DIR}/TB-19-15735-CO394-v2.dna+protein.h5"
  "${H5DIR}/TB22_17721_CO758.dna+protein.h5"
)

# =========================
# Python scripts
# Put these two revised scripts in ${SCRIPTDIR}
# =========================
SCITE_CONVERT_SCRIPT="${SCRIPTDIR}/tapestri_h5_to_scite_dp10_GQ30_AF20_cellQC.py"
VCF_CONVERT_SCRIPT="${SCRIPTDIR}/h5_to_vcf_for_vep_with_dp_vaf_GQ30_AF20_cellQC.py"

if [ ! -f "${SCITE_CONVERT_SCRIPT}" ]; then
    echo "[ERROR] Missing script: ${SCITE_CONVERT_SCRIPT}"
    exit 1
fi

if [ ! -f "${VCF_CONVERT_SCRIPT}" ]; then
    echo "[ERROR] Missing script: ${VCF_CONVERT_SCRIPT}"
    exit 1
fi

# Activate venv if available
if [ -f "${BASE}/scite_venv/bin/activate" ]; then
    source "${BASE}/scite_venv/bin/activate"
fi

echo "=========================================="
echo "[INFO] H5DIR: ${H5DIR}"
echo "[INFO] OUTBASE: ${OUTBASE}"
echo "[INFO] SCITE_OUT: ${SCITE_OUT}"
echo "[INFO] VCF_OUT: ${VCF_OUT}"
echo "[INFO] LOG_OUT: ${LOG_OUT}"
echo "[INFO] MIN_CELL_DEPTH: ${MIN_CELL_DEPTH}"
echo "[INFO] MIN_GENOTYPE_QUALITY: ${MIN_GENOTYPE_QUALITY}"
echo "[INFO] MIN_MUT_AF: ${MIN_MUT_AF}"
echo "[INFO] MAX_MISSING_FRAC: ${MAX_MISSING_FRAC}"
echo "[INFO] MAX_CELL_MISSING_FRAC: ${MAX_CELL_MISSING_FRAC}"
echo "[INFO] MIN_MUTATED_CELLS: ${MIN_MUTATED_CELLS}  # not applied"
echo "=========================================="

for H5 in "${H5_FILES[@]}"; do

    if [ ! -f "${H5}" ]; then
        echo "[WARNING] H5 not found, skipping: ${H5}"
        continue
    fi

    BASENAME=$(basename "${H5}")
    SAMPLE="${BASENAME%.h5}"
    SAMPLE="${SAMPLE%.dna+protein}"

    echo ""
    echo "=========================================="
    echo "[INFO] Processing sample: ${SAMPLE}"
    echo "[INFO] H5: ${H5}"
    echo "=========================================="

    # Check h5 readability and layers
    python - <<PY
import h5py
h5file = "${H5}"
try:
    with h5py.File(h5file, "r") as f:
        print("[OK] H5 opened successfully")
        if "assays/dna_variants/metadata/genome_version" in f:
            gv = f["assays/dna_variants/metadata/genome_version"][()]
            print("[INFO] genome_version:", gv)
        print("[INFO] dna_variants layers:")
        for k in f["assays/dna_variants/layers"].keys():
            d = f["assays/dna_variants/layers"][k]
            print("  ", k, d.shape, d.dtype)
except Exception as e:
    print("[ERROR] Could not open H5:", e)
    raise SystemExit(1)
PY

    if [ $? -ne 0 ]; then
        echo "[SKIP] Broken or incomplete h5: ${H5}"
        continue
    fi

    # ------------------------------------------
    # 1. h5 -> SCITE matrix / geneNames / summary
    # QC rule:
    #   DP < 10 -> missing
    #   GQ < 30 -> missing
    #   NGT 1/2 and 0 < AF < 20 -> missing
    #   NGT 1/2 and AF >= 20 -> mutation positive
    #   AF == 0 -> wild-type if DP/GQ pass
    #   cells with <50% genotypes present -> removed
    #   variants genotyped in <50% retained cells -> removed
    #   MIN_MUTATED_CELLS=5 -> not applied here
    # ------------------------------------------
    echo "[STEP 1] Creating cellQC DP10_GQ30_AF20-filtered SCITE matrix and mutation list"

    python "${SCITE_CONVERT_SCRIPT}" \
        --h5 "${H5}" \
        --outdir "${SCITE_OUT}" \
        --sample "${SAMPLE}" \
        --binary \
        --min-cell-depth "${MIN_CELL_DEPTH}" \
        --min-genotype-quality "${MIN_GENOTYPE_QUALITY}" \
        --min-mut-af "${MIN_MUT_AF}" \
        --max-cell-missing-frac "${MAX_CELL_MISSING_FRAC}" \
        --max-missing-frac "${MAX_MISSING_FRAC}" \
        --min-mutated-cells "${MIN_MUTATED_CELLS}" \
        > "${LOG_OUT}/${SAMPLE}.tapestri_h5_to_scite_DP10_GQ30_AF20_cellQC_noMinMut.log" \
        2> "${LOG_OUT}/${SAMPLE}.tapestri_h5_to_scite_DP10_GQ30_AF20_cellQC_noMinMut.err"

    if [ $? -ne 0 ]; then
        echo "[ERROR] tapestri_h5_to_scite_dp10_GQ30_AF20_cellQC.py failed for ${SAMPLE}"
        cat "${LOG_OUT}/${SAMPLE}.tapestri_h5_to_scite_DP10_GQ30_AF20_cellQC_noMinMut.err"
        continue
    fi

    echo "[DONE] CellQC DP10_GQ30_AF20-filtered SCITE input files created"

    # ------------------------------------------
    # 2. h5 -> all valid REF/ALT VCF for VEP
    #    Summary INFO values are recalculated after cell-level QC.
    # ------------------------------------------
    echo "[STEP 2] Creating all-variant VCF with DP/GQ/AF/cellQC summary for VEP"

    ALL_VCF="${VCF_OUT}/${SAMPLE}.all_variants_for_vep.DP10_GQ30_AF20.cellQC.vcf"
    VCF_CELL_QC="${VCF_OUT}/${SAMPLE}.vcf_cell_qc_summary.tsv"

    python "${VCF_CONVERT_SCRIPT}" \
        --h5 "${H5}" \
        --out "${ALL_VCF}" \
        --min-cell-depth "${MIN_CELL_DEPTH}" \
        --min-genotype-quality "${MIN_GENOTYPE_QUALITY}" \
        --min-mut-af "${MIN_MUT_AF}" \
        --max-cell-missing-frac "${MAX_CELL_MISSING_FRAC}" \
        --cell-qc-out "${VCF_CELL_QC}" \
        > "${LOG_OUT}/${SAMPLE}.h5_to_vcf_for_vep_DP10_GQ30_AF20_cellQC.log" \
        2> "${LOG_OUT}/${SAMPLE}.h5_to_vcf_for_vep_DP10_GQ30_AF20_cellQC.err"

    if [ $? -ne 0 ]; then
        echo "[ERROR] h5_to_vcf_for_vep_with_dp_vaf_GQ30_AF20_cellQC.py failed for ${SAMPLE}"
        cat "${LOG_OUT}/${SAMPLE}.h5_to_vcf_for_vep_DP10_GQ30_AF20_cellQC.err"
        continue
    fi

    echo "[DONE] All-variant cellQC VCF created: ${ALL_VCF}"

    # ------------------------------------------
    # 3. Make matrix-filtered VCF for VEP
    #    Keep only variants retained in mutation_summary.tsv.
    #    This applies cell-level QC + variant missing_fraction <= 0.50.
    #    It does NOT apply MIN_MUTATED_CELLS=5.
    # ------------------------------------------
    echo "[STEP 3] Creating matrix-filtered VCF for VEP"

    SUMMARY="${SCITE_OUT}/${SAMPLE}.mutation_summary.tsv"
    FILTERED_VCF="${VCF_OUT}/${SAMPLE}.matrix_filtered_for_vep.DP10_GQ30_AF20.cellQC.noMinMut.vcf"

    python - <<PY
import pandas as pd
import re

summary_file = "${SUMMARY}"
all_vcf = "${ALL_VCF}"
out_vcf = "${FILTERED_VCF}"

def clean_vcf_id(x):
    x = str(x)
    x = x.replace(" ", "_")
    x = x.replace(";", "_")
    x = x.replace("\\t", "_")
    x = re.sub(r"_+", "_", x)
    return x.strip("_")

summary = pd.read_csv(summary_file, sep="\t", dtype=str)
keep = set(summary["display_name"].map(clean_vcf_id).astype(str))

n_in = 0
n_out = 0

with open(all_vcf) as f, open(out_vcf, "w") as o:
    for line in f:
        if line.startswith("#"):
            o.write(line)
            continue
        n_in += 1
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        vid = parts[2]
        if vid in keep:
            o.write(line)
            n_out += 1

print("[DONE] Wrote:", out_vcf)
print("Input VCF variants:", n_in)
print("Matrix-filtered variants written:", n_out)
print("Mutation summary variants:", len(keep))

if n_out != len(keep):
    print("[WARNING] VCF variants and mutation_summary variants do not fully match.")
    print("[WARNING] Please check whether VCF ID and display_name use the same naming.")
PY

    if [ $? -ne 0 ]; then
        echo "[ERROR] Failed to create matrix-filtered VCF for ${SAMPLE}"
        continue
    fi

    echo "[DONE] Matrix-filtered VCF created: ${FILTERED_VCF}"

    # ------------------------------------------
    # 3b. Optional final VCF filter for the method text:
    #     variants mutated in fewer than 1.0% of cells are removed.
    #     This does NOT affect the SCITE matrix.
    # ------------------------------------------
    FINAL_MUT1_VCF="${VCF_OUT}/${SAMPLE}.final_QC_mutatedPct1.DP10_GQ30_AF20.cellQC.vcf"

    python - <<PY
# Optional final VCF creation.
# Set CREATE_FINAL_MUT1PCT_VCF = True if you want to automatically make this file.
CREATE_FINAL_MUT1PCT_VCF = False
in_vcf = "${FILTERED_VCF}"
out_vcf = "${FINAL_MUT1_VCF}"
min_mut_pct = float("${MIN_MUTATED_PERCENT_FOR_FINAL_VCF}")

import re

def parse_info(info):
    d = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

if CREATE_FINAL_MUT1PCT_VCF:
    n_in = 0
    n_out = 0
    with open(in_vcf) as f, open(out_vcf, "w") as o:
        for line in f:
            if line.startswith("#"):
                o.write(line)
                continue
            n_in += 1
            parts = line.rstrip("\n").split("\t")
            info = parse_info(parts[7])
            mut_pct = float(info.get("MUTATED_PCT_ALL_CELLS_DP10_GQ30_AF20", "0"))
            if mut_pct >= min_mut_pct:
                o.write(line)
                n_out += 1
    print("[DONE] Final mutated >=1% VCF:", out_vcf)
    print("Input variants:", n_in)
    print("Output variants:", n_out)
else:
    print("[SKIP] Optional final mutated >=1% VCF was not created. To create it, set CREATE_FINAL_MUT1PCT_VCF = True in this block.")
PY

    # ------------------------------------------
    # 4. Check output files
    # ------------------------------------------
    echo "[STEP 4] Checking output files"

    echo "[INFO] Retained cells after cell QC:"
    awk -F'\t' 'NR>1 && $6=="True" {n++} END{print n+0}' "${SCITE_OUT}/${SAMPLE}.cell_qc_summary.tsv"

    echo "[INFO] Removed cells by cell QC:"
    awk -F'\t' 'NR>1 && $6!="True" {n++} END{print n+0}' "${SCITE_OUT}/${SAMPLE}.cell_qc_summary.tsv"

    echo "[INFO] SCITE matrix rows = mutations:"
    wc -l "${SCITE_OUT}/${SAMPLE}.scite_matrix.txt"

    echo "[INFO] SCITE matrix columns = cells:"
    awk '{print NF; exit}' "${SCITE_OUT}/${SAMPLE}.scite_matrix.txt"

    echo "[INFO] Matrix unique values:"
    awk '{for(i=1;i<=NF;i++) print $i}' "${SCITE_OUT}/${SAMPLE}.scite_matrix.txt" | sort | uniq -c

    echo "[INFO] geneNames:"
    wc -l "${SCITE_OUT}/${SAMPLE}.geneNames"

    echo "[INFO] mutation_summary:"
    head -n 5 "${SCITE_OUT}/${SAMPLE}.mutation_summary.tsv"

    echo "[INFO] VCF variant counts:"
    echo -n "All valid cellQC VCF variants: "
    grep -v "^#" "${ALL_VCF}" | wc -l

    echo -n "Matrix-filtered VCF variants: "
    grep -v "^#" "${FILTERED_VCF}" | wc -l

    echo "[INFO] Example VCF lines with DP/GQ/AF/cellQC INFO:"
    grep -v "^#" "${FILTERED_VCF}" | head -n 3

    echo "=========================================="
    echo "[DONE] ${SAMPLE}"
    echo "SCITE matrix: ${SCITE_OUT}/${SAMPLE}.scite_matrix.txt"
    echo "geneNames: ${SCITE_OUT}/${SAMPLE}.geneNames"
    echo "mutation summary: ${SCITE_OUT}/${SAMPLE}.mutation_summary.tsv"
    echo "barcodes: ${SCITE_OUT}/${SAMPLE}.barcodes.txt"
    echo "cell QC summary: ${SCITE_OUT}/${SAMPLE}.cell_qc_summary.tsv"
    echo "clone pattern summary: ${SCITE_OUT}/${SAMPLE}.clone_pattern_summary.tsv"
    echo "all VCF: ${ALL_VCF}"
    echo "matrix-filtered VCF: ${FILTERED_VCF}"
    echo "=========================================="

done

echo ""
echo "=========================================="
echo "[ALL DONE]"
echo "Output directory:"
echo "${OUTBASE}"
echo "=========================================="
