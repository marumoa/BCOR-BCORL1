#!/usr/bin/env python3
"""
Create SCITE input files from Mission Bio / Tapestri h5 using the QC rule:
  DP < 10                         -> missing
  GQ < 30                         -> missing
  NGT in 1,2 and 0 < AF < 20%      -> missing
  NGT in 1,2 and AF >= 20%         -> mutation positive
  AF == 0 with sufficient DP/GQ     -> wild-type / mutation negative
  NGT == 0 with sufficient DP/GQ    -> wild-type / mutation negative

Additional QC implemented here:
  cell-level QC: remove cells with <50% genotypes present
                 i.e. cell missing_fraction > 0.50
  variant-level QC: remove variants genotyped in <50% retained cells
                    i.e. variant missing_fraction > 0.50

By default, MIN_MUTATED_CELLS is NOT applied. If you want it later, set
--min-mutated-cells to a positive integer.
"""

import argparse
import re
from collections import Counter
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


def safe_mean(x):
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return np.nan
    return float(np.nanmean(x))


def safe_median(x):
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return np.nan
    return float(np.nanmedian(x))


def fmt_num(x, digits=4):
    try:
        x = float(x)
    except Exception:
        return ""
    if np.isnan(x) or np.isinf(x):
        return ""
    return f"{x:.{digits}f}"


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

    return {
        "high_quality": high_quality,
        "mutation_positive": mutation_positive,
        "wildtype": wildtype,
        "called": called,
        "final_missing": final_missing,
        "low_af_missing": low_af_missing,
        "af_zero_wt": af_zero_wt,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--binary", action="store_true", help="Kept for compatibility. Output is 0/1/3 SCITE matrix.")
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
    parser.add_argument("--max-missing-frac", type=float, default=0.50,
                        help="Remove variants with missing fraction above this threshold after cell QC.")
    parser.add_argument("--min-mutated-cells", type=int, default=0,
                        help="Default 0 = not applied. Use positive integer if needed.")
    parser.add_argument("--min-mutated-percent", type=float, default=0.0,
                        help="Default 0 = not applied. Use 1.0 to apply variants mutated in >=1%% retained cells.")
    parser.add_argument("--missing-code", default="3", help="Missing genotype code in SCITE matrix. Default: 3.")
    parser.add_argument("--keep-chr", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mutated_genotypes = parse_int_list(args.mutated_genotypes)
    wildtype_genotypes = parse_int_list(args.wildtype_genotypes)

    with h5py.File(args.h5, "r") as f:
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
    if args.passing_only and filtered is not None:
        eligible_variant &= make_filter_pass_mask(filtered)

    # Cell-level QC: use eligible variants only. Present genotype = final called non-missing.
    called_for_cell_qc = calls["called"][eligible_variant, :]
    if called_for_cell_qc.shape[0] == 0:
        raise RuntimeError("No eligible variants available for cell QC.")

    cell_called = called_for_cell_qc.sum(axis=0).astype(int)
    cell_total = int(called_for_cell_qc.shape[0])
    cell_present_fraction = cell_called / cell_total
    cell_missing_fraction = 1.0 - cell_present_fraction
    keep_cell = cell_missing_fraction <= args.max_cell_missing_frac

    # Variant-level QC after cell QC.
    called_ret = calls["called"][:, keep_cell]
    mut_ret = calls["mutation_positive"][:, keep_cell]
    wt_ret = calls["wildtype"][:, keep_cell]
    miss_ret = calls["final_missing"][:, keep_cell]

    retained_cell_count = int(keep_cell.sum())
    if retained_cell_count == 0:
        raise RuntimeError("No cells retained after cell-level missing fraction filter.")

    called_cells = called_ret.sum(axis=1).astype(int)
    mutated_cells = mut_ret.sum(axis=1).astype(int)
    wildtype_cells = wt_ret.sum(axis=1).astype(int)
    missing_cells = miss_ret.sum(axis=1).astype(int)
    missing_fraction = missing_cells / retained_cell_count
    mutated_percent = mutated_cells / retained_cell_count * 100.0
    mutated_percent_called = np.divide(
        mutated_cells, called_cells,
        out=np.zeros_like(mutated_cells, dtype=float),
        where=called_cells > 0,
    ) * 100.0

    keep_variant = eligible_variant & (missing_fraction <= args.max_missing_frac)
    if args.min_mutated_cells and args.min_mutated_cells > 0:
        keep_variant &= mutated_cells >= args.min_mutated_cells
    if args.min_mutated_percent and args.min_mutated_percent > 0:
        keep_variant &= mutated_percent >= args.min_mutated_percent

    keep_variant_idx = np.where(keep_variant)[0]

    # Matrix: rows = retained variants, columns = retained cells.
    matrix = np.full((len(keep_variant_idx), retained_cell_count), args.missing_code, dtype=object)
    for row_i, var_i in enumerate(keep_variant_idx):
        matrix[row_i, mut_ret[var_i, :]] = "1"
        matrix[row_i, wt_ret[var_i, :]] = "0"

    matrix_file = outdir / f"{args.sample}.scite_matrix.txt"
    gene_file = outdir / f"{args.sample}.geneNames"
    summary_file = outdir / f"{args.sample}.mutation_summary.tsv"
    barcodes_file = outdir / f"{args.sample}.barcodes.txt"
    cell_qc_file = outdir / f"{args.sample}.cell_qc_summary.tsv"
    clone_file = outdir / f"{args.sample}.clone_pattern_summary.tsv"

    with open(matrix_file, "w") as f:
        for row in matrix:
            f.write(" ".join(map(str, row)) + "\n")

    with open(gene_file, "w") as f:
        for i in keep_variant_idx:
            f.write(clean_vcf_id(vid[i]) + "\n")

    with open(barcodes_file, "w") as f:
        for bc, keep in zip(barcodes, keep_cell):
            if keep:
                f.write(str(bc) + "\n")

    cell_df = pd.DataFrame({
        "barcode": barcodes,
        "total_variants_for_cell_qc": cell_total,
        "called_genotypes": cell_called,
        "genotype_present_fraction": cell_present_fraction,
        "missing_fraction": cell_missing_fraction,
        "keep_cell": keep_cell,
    })
    cell_df.to_csv(cell_qc_file, sep="\t", index=False)

    summary_rows = []
    for i in keep_variant_idx:
        called_mask = called_ret[i, :]
        mut_mask = mut_ret[i, :]
        summary_rows.append({
            "display_name": vid[i],
            "vcf_id": clean_vcf_id(vid[i]),
            "chrom": chrom[i],
            "pos": int(pos[i]),
            "ref": ref[i],
            "alt": alt[i],
            "total_cells_original": n_cell,
            "total_cells_after_cell_qc": retained_cell_count,
            "removed_cells_by_cell_qc": int(n_cell - retained_cell_count),
            "cell_qc_max_missing_fraction": args.max_cell_missing_frac,
            "called_cells": int(called_cells[i]),
            "wildtype_cells": int(wildtype_cells[i]),
            "mutated_cells": int(mutated_cells[i]),
            "mutated_percent": fmt_num(mutated_percent[i]),
            "mutated_percent_called_cells": fmt_num(mutated_percent_called[i]),
            "missing_cells": int(missing_cells[i]),
            "missing_fraction": fmt_num(missing_fraction[i]),
            "mean_dp_called": fmt_num(safe_mean(dp[i, keep_cell][called_mask])),
            "median_dp_called": fmt_num(safe_median(dp[i, keep_cell][called_mask])),
            "mean_gq_called": fmt_num(safe_mean(gq[i, keep_cell][called_mask])),
            "median_gq_called": fmt_num(safe_median(gq[i, keep_cell][called_mask])),
            "mean_af_mut": fmt_num(safe_mean(af[i, keep_cell][mut_mask])),
            "median_af_mut": fmt_num(safe_median(af[i, keep_cell][mut_mask])),
            "min_dp_cutoff": args.min_cell_depth,
            "min_gq_cutoff": args.min_genotype_quality,
            "min_af_mut_cutoff_percent": args.min_mut_af,
            "af_scale_used": detected_af_scale,
            "af_cutoff_in_h5_scale": af_cutoff,
            "variant_max_missing_fraction_filter": args.max_missing_frac,
            "min_mutated_cells_filter": args.min_mutated_cells,
            "min_mutated_percent_filter": args.min_mutated_percent,
        })
    pd.DataFrame(summary_rows).to_csv(summary_file, sep="\t", index=False)

    # Cell clone/pattern summary across retained variants.
    patterns = []
    retained_barcodes = [bc for bc, keep in zip(barcodes, keep_cell) if keep]
    if matrix.shape[0] > 0 and matrix.shape[1] > 0:
        matrix_t = matrix.T
        for bc, vals in zip(retained_barcodes, matrix_t):
            patterns.append(("".join(map(str, vals)), bc))
    cnt = Counter(p for p, _ in patterns)
    with open(clone_file, "w") as f:
        f.write("pattern\tcell_count\n")
        for pattern, count in cnt.most_common():
            f.write(f"{pattern}\t{count}\n")

    print(f"[DONE] Wrote SCITE matrix: {matrix_file}")
    print(f"[DONE] Wrote geneNames: {gene_file}")
    print(f"[DONE] Wrote mutation summary: {summary_file}")
    print(f"[DONE] Wrote retained barcodes: {barcodes_file}")
    print(f"[DONE] Wrote cell QC summary: {cell_qc_file}")
    print(f"Original cells: {n_cell}")
    print(f"Retained cells after cell-level QC: {retained_cell_count}")
    print(f"Removed cells by cell-level QC: {n_cell - retained_cell_count}")
    print(f"Eligible variants for cell QC: {int(eligible_variant.sum())}")
    print(f"Retained variants after variant-level QC: {len(keep_variant_idx)}")
    print(f"MIN_MUTATED_CELLS applied: {args.min_mutated_cells}")
    print(f"MIN_MUTATED_PERCENT applied: {args.min_mutated_percent}")


if __name__ == "__main__":
    main()
