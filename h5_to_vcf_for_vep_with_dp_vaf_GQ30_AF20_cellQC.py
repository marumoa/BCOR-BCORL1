#!/usr/bin/env python3
"""
Convert Mission Bio / Tapestri h5 DNA variants to a VCF for VEP annotation.

This version adds cell-level QC:
  remove cells with <50% genotypes present after DP/GQ/AF rule
  i.e. remove cells with missing_fraction > --max-cell-missing-frac.

The VCF still writes one record per eligible h5 variant with valid CHROM/POS/REF/ALT.
Cell-level summary INFO values are recalculated after removing low-quality cells.
"""

import argparse
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def dec(x):
    if isinstance(x, (bytes, bytearray)):
        return x.decode("utf-8", errors="replace")
    return str(x)


def decode_array(x):
    x = np.asarray(x)
    if x.dtype.kind in {"S", "O"}:
        return np.array([
            v.decode("utf-8", errors="replace") if isinstance(v, (bytes, bytearray)) else str(v)
            for v in x
        ])
    return x


def clean_vcf_id(x):
    x = str(x)
    x = x.replace(" ", "_").replace(";", "_").replace("\t", "_")
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def first_attr(attr_group, names):
    if attr_group is None:
        return None
    for n in names:
        if n in attr_group:
            return decode_array(attr_group[n][()])
    return None


def find_layer(h5, requested_layer="NGT"):
    candidates = []
    target_layers = {"NGT", "NGT_FILTERED"} if requested_layer == "auto" else {requested_layer}

    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset):
            base = name.split("/")[-1]
            if base in target_layers and "/layers/" in name:
                candidates.append(name)

    h5.visititems(visitor)
    if not candidates:
        raise RuntimeError(f"{requested_layer} layer was not found under any /layers/ path.")

    def score(path):
        s = 0
        lower = path.lower()
        if "dna" in lower:
            s += 10
        if path.endswith("/NGT"):
            s += 5
        if path.endswith("/NGT_FILTERED"):
            s += 3
        return s

    return sorted(candidates, key=score, reverse=True)[0]


def load_assay_attrs(h5, assay_path):
    assay = h5[assay_path]
    ca = None
    ra = None
    for c in ["ca", "col_attrs", "column_attrs"]:
        if c in assay:
            ca = assay[c]
            break
    for r in ["ra", "row_attrs", "row_attributes"]:
        if r in assay:
            ra = assay[r]
            break
    return ca, ra


def layer_if_exists(h5, assay_path, layer_name):
    p = assay_path + "/layers/" + layer_name
    if p in h5:
        return np.asarray(h5[p][()])
    return None


def orient_matrix(arr, n_cell, n_var):
    if arr is None:
        return None
    if arr.shape == (n_cell, n_var):
        return arr.T
    if arr.shape == (n_var, n_cell):
        return arr
    raise RuntimeError(f"Layer orientation is unclear. layer={arr.shape}, cells={n_cell}, variants={n_var}")


def vector_for_variant(arr, i):
    if arr is None:
        return None
    return arr[i, :]


def safe_float(x):
    if x is None:
        return None
    try:
        x = float(x)
        if np.isnan(x) or np.isinf(x):
            return None
        return x
    except Exception:
        return None


def fmt_float(x, ndigit=4):
    x = safe_float(x)
    if x is None:
        return None
    return f"{x:.{ndigit}f}"


def parse_int_list(x):
    out = []
    for v in str(x).split(","):
        v = v.strip()
        if v:
            out.append(int(v))
    return out


def make_filter_pass_mask(filter_values):
    arr = np.asarray(filter_values)
    if arr.dtype.kind in {"i", "u", "f", "b"}:
        return arr == 0
    arr = decode_array(arr)
    arr_str = np.array([str(x).strip().lower() for x in arr])
    return np.isin(arr_str, ["0", "false", "pass", "passed", "none", "nan", ""])


def infer_af_cutoff(af_float, min_mut_af_percent, af_scale):
    valid = af_float[np.isfinite(af_float)]
    valid = valid[valid >= 0]
    if af_scale == "fraction":
        return min_mut_af_percent / 100.0, "fraction"
    if af_scale == "percent":
        return min_mut_af_percent, "percent"
    if valid.size == 0:
        return min_mut_af_percent, "percent"
    vmax = float(np.nanmax(valid))
    if vmax <= 1.5:
        return min_mut_af_percent / 100.0, "fraction"
    return min_mut_af_percent, "percent"


def strip_chr_prefix(chrom):
    c = str(chrom)
    return c[3:] if c.lower().startswith("chr") else c


def is_valid_ref_alt(ref, alt):
    ref = str(ref)
    alt = str(alt)
    if ref in ["", ".", "nan", "None"]:
        return False
    if alt in ["", ".", "nan", "None"]:
        return False
    return True


def build_calls(ngt, dp, gq, af, min_cell_depth, min_genotype_quality, af_cutoff,
                mutated_genotypes, wildtype_genotypes):
    high_quality = (
        np.isfinite(dp)
        & np.isfinite(gq)
        & (dp >= min_cell_depth)
        & (gq >= min_genotype_quality)
    )
    af_valid = np.isfinite(af) & (af >= 0)
    mutated_ngt = np.isin(ngt, mutated_genotypes)
    wildtype_ngt = np.isin(ngt, wildtype_genotypes)

    mutation_positive = high_quality & mutated_ngt & af_valid & (af >= af_cutoff)
    low_af_missing = high_quality & mutated_ngt & af_valid & (af > 0) & (af < af_cutoff)
    af_zero_wt = high_quality & af_valid & (af == 0)
    wildtype = (high_quality & wildtype_ngt) | af_zero_wt
    called = mutation_positive | wildtype
    final_missing = ~called
    low_quality_missing = ~high_quality

    return {
        "high_quality": high_quality,
        "mutation_positive": mutation_positive,
        "wildtype": wildtype,
        "called": called,
        "final_missing": final_missing,
        "low_af_missing": low_af_missing,
        "af_zero_wt": af_zero_wt,
        "low_quality_missing": low_quality_missing,
    }


def fmt_info_float(key, value, ndigit=4):
    val = fmt_float(value, ndigit)
    return f"{key}={val}" if val is not None else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--ngt-layer", choices=["NGT", "NGT_FILTERED", "auto"], default="NGT")
    parser.add_argument("--passing-only", action="store_true")
    parser.add_argument("--min-cell-depth", type=int, default=10)
    parser.add_argument("--min-genotype-quality", type=float, default=30)
    parser.add_argument("--min-mut-af", type=float, default=20.0)
    parser.add_argument("--af-scale", choices=["auto", "percent", "fraction"], default="auto")
    parser.add_argument("--mutated-genotypes", default="1,2")
    parser.add_argument("--wildtype-genotypes", default="0")
    parser.add_argument("--max-cell-missing-frac", type=float, default=0.50,
                        help="Remove cells with missing fraction above this threshold after DP/GQ/AF rule.")
    parser.add_argument("--keep-chr", action="store_true")
    parser.add_argument("--cell-qc-out", default=None,
                        help="Optional TSV path to write cell-level QC metrics.")
    args = parser.parse_args()

    h5_path = Path(args.h5)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    mutated_genotypes = parse_int_list(args.mutated_genotypes)
    wildtype_genotypes = parse_int_list(args.wildtype_genotypes)

    with h5py.File(h5_path, "r") as f:
        ngt_layer_path = find_layer(f, requested_layer=args.ngt_layer)
        assay_path = ngt_layer_path.rsplit("/layers/", 1)[0]
        ca, ra = load_assay_attrs(f, assay_path)
        if ca is None:
            raise RuntimeError(f"Column attributes/ca were not found under assay: {assay_path}")
        if ra is None:
            raise RuntimeError(f"Row attributes/ra were not found under assay: {assay_path}")

        chrom_arr = first_attr(ca, ["CHROM", "chrom", "chr", "Chromosome"])
        pos_arr = first_attr(ca, ["POS", "pos", "position", "Position"])
        ref_arr = first_attr(ca, ["REF", "ref"])
        alt_arr = first_attr(ca, ["ALT", "alt"])
        id_arr = first_attr(ca, ["id", "ids", "ID", "variant_id"])
        barcodes = first_attr(ra, ["barcode", "barcodes", "id", "ids", "cell_barcode"])
        filtered = first_attr(ca, ["filtered", "FILTERED"])

        if chrom_arr is None or pos_arr is None or ref_arr is None or alt_arr is None:
            raise RuntimeError("CHROM/POS/REF/ALT were not found in h5 variant attributes.")
        if id_arr is None:
            raise RuntimeError("Variant IDs were not found in h5 variant attributes.")
        if barcodes is None:
            raise RuntimeError("Cell barcodes were not found in h5 row attributes.")

        n_var = len(id_arr)
        n_cell = len(barcodes)
        chrom = [dec(x) if args.keep_chr else strip_chr_prefix(dec(x)) for x in chrom_arr]
        pos = pos_arr
        ref = [dec(x) for x in ref_arr]
        alt = [dec(x) for x in alt_arr]
        vid = [dec(x) for x in id_arr]
        barcodes = [dec(x) for x in barcodes]

        ngt = orient_matrix(np.asarray(f[ngt_layer_path][()]), n_cell, n_var)
        dp_raw = layer_if_exists(f, assay_path, "DP")
        gq_raw = layer_if_exists(f, assay_path, "GQ")
        af_raw = layer_if_exists(f, assay_path, "AF")
        if dp_raw is None:
            raise RuntimeError("DP layer was not found.")
        if gq_raw is None:
            raise RuntimeError("GQ layer was not found.")
        if af_raw is None:
            raise RuntimeError("AF layer was not found. Cannot apply AF >= 20% rule.")
        dp = orient_matrix(dp_raw, n_cell, n_var).astype(float)
        gq = orient_matrix(gq_raw, n_cell, n_var).astype(float)
        af = orient_matrix(af_raw, n_cell, n_var).astype(float)

    af_cutoff, detected_af_scale = infer_af_cutoff(af, args.min_mut_af, args.af_scale)
    calls = build_calls(
        ngt, dp, gq, af,
        args.min_cell_depth,
        args.min_genotype_quality,
        af_cutoff,
        mutated_genotypes,
        wildtype_genotypes,
    )

    valid_ref_alt = np.array([is_valid_ref_alt(r, a) for r, a in zip(ref, alt)], dtype=bool)
    eligible_variant = valid_ref_alt.copy()
    variant_pass_mask = None
    if filtered is not None:
        variant_pass_mask = make_filter_pass_mask(filtered)
        if args.passing_only:
            eligible_variant &= variant_pass_mask

    # Cell-level QC: use eligible variants only.
    called_for_cell_qc = calls["called"][eligible_variant, :]
    if called_for_cell_qc.shape[0] == 0:
        raise RuntimeError("No eligible variants available for cell QC.")

    cell_called = called_for_cell_qc.sum(axis=0).astype(int)
    cell_total = int(called_for_cell_qc.shape[0])
    cell_present_fraction = cell_called / cell_total
    cell_missing_fraction = 1.0 - cell_present_fraction
    keep_cell = cell_missing_fraction <= args.max_cell_missing_frac

    if args.cell_qc_out:
        cell_qc_path = Path(args.cell_qc_out)
    else:
        cell_qc_path = out_path.with_suffix(".cell_qc_summary.tsv")
    cell_qc_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "barcode": barcodes,
        "total_variants_for_cell_qc": cell_total,
        "called_genotypes": cell_called,
        "genotype_present_fraction": cell_present_fraction,
        "missing_fraction": cell_missing_fraction,
        "keep_cell": keep_cell,
    }).to_csv(cell_qc_path, sep="\t", index=False)

    n_written = 0
    n_skipped_filtered = 0
    n_skipped_invalid_ref_alt = 0
    retained_cell_count = int(keep_cell.sum())
    if retained_cell_count == 0:
        raise RuntimeError("No cells retained after cell-level missing fraction filter.")

    with open(out_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=Tapestri_h5_to_VCF_for_VEP_DP10_GQ30_AF20_cellQC\n")
        out.write(f'##INFO=<ID=H5_NGT_LAYER,Number=1,Type=String,Description="NGT layer used: {args.ngt_layer}">\n')
        out.write('##INFO=<ID=H5_FILTERED,Number=1,Type=Integer,Description="Variant-level filtered flag from h5 ca/filtered; 0=pass, 1=filtered if numeric">\n')
        out.write('##INFO=<ID=TOTAL_CELLS_ORIGINAL,Number=1,Type=Integer,Description="Original total number of cells in h5">\n')
        out.write('##INFO=<ID=TOTAL_CELLS_AFTER_CELL_QC,Number=1,Type=Integer,Description="Cells retained after removing cells with excessive missing genotypes">\n')
        out.write('##INFO=<ID=REMOVED_CELLS_BY_CELL_QC,Number=1,Type=Integer,Description="Cells removed because genotype present fraction was less than the threshold">\n')
        out.write('##INFO=<ID=CELL_QC_MAX_MISSING_FRACTION,Number=1,Type=Float,Description="Cell-level missing fraction cutoff; cells above this are removed">\n')
        out.write('##INFO=<ID=HIGH_QUALITY_CELLS_DP10_GQ30,Number=1,Type=Integer,Description="Retained cells with DP >= min-cell-depth and GQ >= min-genotype-quality">\n')
        out.write('##INFO=<ID=CALLED_CELLS_DP10_GQ30_AF20RULE,Number=1,Type=Integer,Description="Retained cells with final non-missing call after DP/GQ/AF rule">\n')
        out.write('##INFO=<ID=MUTATED_CELLS_DP10_GQ30_AF20,Number=1,Type=Integer,Description="Retained cells with NGT in mutated genotypes, DP/GQ pass, and AF >= AF cutoff">\n')
        out.write('##INFO=<ID=MUTATED_PCT_ALL_CELLS_DP10_GQ30_AF20,Number=1,Type=Float,Description="Percent of retained cells that are mutation positive by DP/GQ/AF rule">\n')
        out.write('##INFO=<ID=MUTATED_PCT_CALLED_CELLS_DP10_GQ30_AF20,Number=1,Type=Float,Description="Percent of called non-missing retained cells that are mutation positive by DP/GQ/AF rule">\n')
        out.write('##INFO=<ID=WILDTYPE_CELLS_DP10_GQ30_AF20RULE,Number=1,Type=Integer,Description="Retained cells called wild-type by DP/GQ/AF rule">\n')
        out.write('##INFO=<ID=AF_ZERO_WT_CELLS,Number=1,Type=Integer,Description="Retained cells with DP/GQ pass and AF == 0; called wild-type">\n')
        out.write('##INFO=<ID=LOW_QUALITY_MISSING_CELLS,Number=1,Type=Integer,Description="Retained cells made missing because DP or GQ failed">\n')
        out.write('##INFO=<ID=LOW_AF_MISSING_CELLS,Number=1,Type=Integer,Description="Retained cells with mutated NGT and 0 < AF < AF cutoff; made missing">\n')
        out.write('##INFO=<ID=MISSING_CELLS_FINAL,Number=1,Type=Integer,Description="Final missing retained cells after DP/GQ/AF rule">\n')
        out.write('##INFO=<ID=MISSING_FRACTION_FINAL,Number=1,Type=Float,Description="Final missing fraction among retained cells after DP/GQ/AF rule">\n')
        out.write('##INFO=<ID=MEAN_DP_CALLED,Number=1,Type=Float,Description="Mean DP among final called retained cells">\n')
        out.write('##INFO=<ID=MEDIAN_DP_CALLED,Number=1,Type=Float,Description="Median DP among final called retained cells">\n')
        out.write('##INFO=<ID=MEAN_GQ_CALLED,Number=1,Type=Float,Description="Mean GQ among final called retained cells">\n')
        out.write('##INFO=<ID=MEDIAN_GQ_CALLED,Number=1,Type=Float,Description="Median GQ among final called retained cells">\n')
        out.write('##INFO=<ID=MEAN_AF_MUT,Number=1,Type=Float,Description="Mean AF among mutation-positive retained cells; same scale as h5 AF layer">\n')
        out.write('##INFO=<ID=MEDIAN_AF_MUT,Number=1,Type=Float,Description="Median AF among mutation-positive retained cells; same scale as h5 AF layer">\n')
        out.write('##INFO=<ID=MIN_DP_CUTOFF,Number=1,Type=Integer,Description="DP cutoff used">\n')
        out.write('##INFO=<ID=MIN_GQ_CUTOFF,Number=1,Type=Float,Description="GQ cutoff used">\n')
        out.write('##INFO=<ID=MIN_AF_MUT_CUTOFF_PERCENT,Number=1,Type=Float,Description="AF mutation cutoff requested in percent">\n')
        out.write('##INFO=<ID=AF_SCALE_USED,Number=1,Type=String,Description="AF scale detected or specified: percent or fraction">\n')
        out.write('##INFO=<ID=AF_CUTOFF_IN_H5_SCALE,Number=1,Type=Float,Description="AF cutoff converted to h5 AF layer scale">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for i, (c, p, r, a, v) in enumerate(zip(chrom, pos, ref, alt, vid)):
            if args.passing_only and variant_pass_mask is not None and not bool(variant_pass_mask[i]):
                n_skipped_filtered += 1
                continue
            if not is_valid_ref_alt(r, a):
                n_skipped_invalid_ref_alt += 1
                continue

            high_quality = calls["high_quality"][i, keep_cell]
            mutation_positive = calls["mutation_positive"][i, keep_cell]
            wildtype = calls["wildtype"][i, keep_cell]
            called = calls["called"][i, keep_cell]
            low_af_missing = calls["low_af_missing"][i, keep_cell]
            low_quality_missing = calls["low_quality_missing"][i, keep_cell]
            af_zero_wt = calls["af_zero_wt"][i, keep_cell]
            final_missing = calls["final_missing"][i, keep_cell]

            dp_i = dp[i, keep_cell]
            gq_i = gq[i, keep_cell]
            af_i = af[i, keep_cell]

            high_quality_cells = int(high_quality.sum())
            called_cells = int(called.sum())
            mutated_cells = int(mutation_positive.sum())
            wildtype_cells = int(wildtype.sum())
            low_quality_missing_cells = int(low_quality_missing.sum())
            low_af_missing_cells = int(low_af_missing.sum())
            af_zero_cells = int(af_zero_wt.sum())
            missing_cells_final = int(final_missing.sum())
            missing_fraction_final = missing_cells_final / retained_cell_count if retained_cell_count else 0.0
            mutated_pct_all = mutated_cells / retained_cell_count * 100 if retained_cell_count else 0.0
            mutated_pct_called = mutated_cells / called_cells * 100 if called_cells else 0.0

            info_items = []
            if filtered is not None:
                try:
                    info_items.append(f"H5_FILTERED={int(filtered[i])}")
                except Exception:
                    pass
            info_items.append(f"H5_NGT_LAYER={args.ngt_layer}")
            info_items.append(f"TOTAL_CELLS_ORIGINAL={n_cell}")
            info_items.append(f"TOTAL_CELLS_AFTER_CELL_QC={retained_cell_count}")
            info_items.append(f"REMOVED_CELLS_BY_CELL_QC={n_cell - retained_cell_count}")
            info_items.append(f"CELL_QC_MAX_MISSING_FRACTION={args.max_cell_missing_frac:.4f}")
            info_items.append(f"HIGH_QUALITY_CELLS_DP10_GQ30={high_quality_cells}")
            info_items.append(f"CALLED_CELLS_DP10_GQ30_AF20RULE={called_cells}")
            info_items.append(f"MUTATED_CELLS_DP10_GQ30_AF20={mutated_cells}")
            info_items.append(f"MUTATED_PCT_ALL_CELLS_DP10_GQ30_AF20={mutated_pct_all:.4f}")
            info_items.append(f"MUTATED_PCT_CALLED_CELLS_DP10_GQ30_AF20={mutated_pct_called:.4f}")
            info_items.append(f"WILDTYPE_CELLS_DP10_GQ30_AF20RULE={wildtype_cells}")
            info_items.append(f"AF_ZERO_WT_CELLS={af_zero_cells}")
            info_items.append(f"LOW_QUALITY_MISSING_CELLS={low_quality_missing_cells}")
            info_items.append(f"LOW_AF_MISSING_CELLS={low_af_missing_cells}")
            info_items.append(f"MISSING_CELLS_FINAL={missing_cells_final}")
            info_items.append(f"MISSING_FRACTION_FINAL={missing_fraction_final:.4f}")
            info_items.append(f"MIN_DP_CUTOFF={args.min_cell_depth}")
            info_items.append(f"MIN_GQ_CUTOFF={args.min_genotype_quality}")
            info_items.append(f"MIN_AF_MUT_CUTOFF_PERCENT={args.min_mut_af}")
            info_items.append(f"AF_SCALE_USED={detected_af_scale}")
            info_items.append(f"AF_CUTOFF_IN_H5_SCALE={af_cutoff:.6g}")

            if called_cells > 0:
                for key, val in [
                    ("MEAN_DP_CALLED", np.nanmean(dp_i[called])),
                    ("MEDIAN_DP_CALLED", np.nanmedian(dp_i[called])),
                    ("MEAN_GQ_CALLED", np.nanmean(gq_i[called])),
                    ("MEDIAN_GQ_CALLED", np.nanmedian(gq_i[called])),
                ]:
                    item = fmt_info_float(key, val, 4)
                    if item:
                        info_items.append(item)
            if mutated_cells > 0:
                af_mut = af_i[mutation_positive].astype(float)
                for key, val in [
                    ("MEAN_AF_MUT", np.nanmean(af_mut)),
                    ("MEDIAN_AF_MUT", np.nanmedian(af_mut)),
                ]:
                    item = fmt_info_float(key, val, 4)
                    if item:
                        info_items.append(item)

            info = ";".join(info_items) if info_items else "."
            out.write(f"{c}\t{int(p)}\t{clean_vcf_id(v)}\t{r}\t{a}\t.\tPASS\t{info}\n")
            n_written += 1

    print(f"[DONE] Wrote VCF: {out_path}")
    print(f"[DONE] Wrote cell QC summary: {cell_qc_path}")
    print(f"Original cells: {n_cell}")
    print(f"Retained cells after cell-level QC: {retained_cell_count}")
    print(f"Removed cells by cell-level QC: {n_cell - retained_cell_count}")
    print(f"Eligible variants for cell QC: {int(eligible_variant.sum())}")
    print(f"Variants written: {n_written}")
    print(f"Skipped by --passing-only: {n_skipped_filtered}")
    print(f"Skipped due to invalid REF/ALT: {n_skipped_invalid_ref_alt}")
    print(f"NGT layer requested: {args.ngt_layer}")
    print(f"DP cutoff: >= {args.min_cell_depth}")
    print(f"GQ cutoff: >= {args.min_genotype_quality}")
    print(f"AF mutation cutoff requested: >= {args.min_mut_af}%")
    print(f"AF scale used: {detected_af_scale}")
    print(f"AF cutoff in h5 scale: >= {af_cutoff}")
    print(f"Cell-level max missing fraction: <= {args.max_cell_missing_frac}")


if __name__ == "__main__":
    main()
