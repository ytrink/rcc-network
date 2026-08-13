#!/usr/bin/env python3
"""
utils/scenicplus_to_bedpe.py

Converts SCENIC+ eRegulon output to BEDPE format for pyGenomeTracks arc visualization.
Also optionally writes a BED file of enhancer regions.

Usage:
    python scenicplus_to_bedpe.py \
        --ereg data/eRegulon_direct.tsv \
        --tss data/tss.bed \
        --out output/links_SPI1_TREM2.bedpe \
        --out-enhancers-bed output/enh_SPI1_TREM2.bed \
        --tf SPI1 --gene TREM2 \
        --region chr6:40963102-41363102 \
        --top-n 200 --cmap Blues --log-score --merge-enhancers
"""

import re
import argparse

import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# ── Constants ──────────────────────────────────────────────────────────────────
REGION_REGEX = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")

# ── Helpers ────────────────────────────────────────────────────────────────────
def parse_region(s: str):
    """Parse a genomic region string (e.g. chr6:1000-2000) into (chr, start, end)."""
    m = REGION_REGEX.match(str(s))
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def score_to_rgb(scores: pd.Series, cmap_name: str, vmin=None, vmax=None):
    """
    Map numeric scores to RGB strings 'r,g,b' (0-255) using a matplotlib colormap.

    Returns:
        rgb_series : pd.Series of 'r,g,b' strings
        vmin, vmax : floats — the color scale bounds used
    """
    cmap = cm.get_cmap(cmap_name)
    vals = scores.to_numpy(dtype=float)

    if vmin is None:
        vmin = np.nanmin(vals) if len(vals) else 0.0
    if vmax is None:
        vmax = np.nanmax(vals) if len(vals) else 1.0
    if vmax == vmin:
        vmax = vmin + 1e-9

    norm  = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    rgba  = cmap(norm(vals))
    rgb255 = (rgba[:, :3] * 255).astype(int)

    rgb_series = pd.Series(
        [f"{r},{g},{b}" for r, g, b in rgb255],
        index=scores.index
    )
    return rgb_series, vmin, vmax


def load_tss(path: str) -> pd.DataFrame:
    """
    Load a TSS BED file (>=4 cols: chr start end gene).
    Returns a DataFrame indexed by gene name with chr2, start2, end2 columns.
    """
    tss = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    if tss.shape[1] < 4:
        raise ValueError("tss.bed must have at least 4 columns: chr start end gene")

    tss = tss.iloc[:, :4].copy()
    tss.columns = ["chr2", "start_raw", "end_raw", "Gene"]
    tss["tss"]    = ((tss["start_raw"].astype(int) + tss["end_raw"].astype(int)) // 2).astype(int)
    tss["start2"] = tss["tss"]
    tss["end2"]   = tss["tss"] + 1

    return tss.set_index("Gene")[["chr2", "start2", "end2"]]


def load_eregulon(path: str, tf: str = None, gene: str = None) -> pd.DataFrame:
    """
    Load SCENIC+ eRegulon TSV and optionally filter by TF and/or gene.
    """
    df = pd.read_csv(path, sep="\t")

    for col in ["Region", "Gene"]:
        if col not in df.columns:
            raise ValueError(f"Missing required column in eRegulon file: {col}")

    if gene is not None:
        df = df[df["Gene"] == gene]
    if tf is not None and "TF" in df.columns:
        df = df[df["TF"] == tf]

    return df


def filter_by_locus(df: pd.DataFrame, region: str) -> pd.DataFrame:
    """Keep only rows where enhancer or TSS overlaps the given locus."""
    m = REGION_REGEX.match(region)
    if not m:
        raise ValueError("--region must be formatted as chr12:49200000-49300000")

    locus_chr   = m.group(1)
    locus_start = int(m.group(2))
    locus_end   = int(m.group(3))

    in_enh = (
        (df["chr1"] == locus_chr) &
        (df["end1"] >= locus_start) &
        (df["start1"] <= locus_end)
    )
    in_tss = (
        (df["chr2"] == locus_chr) &
        (df["end2"] >= locus_start) &
        (df["start2"] <= locus_end)
    )
    return df[in_enh | in_tss].copy()


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Convert SCENIC+ eRegulon output to BEDPE for pyGenomeTracks."
    )
    ap.add_argument("--ereg",             required=True,  help="eRegulon_direct.tsv or eRegulon_extended.tsv")
    ap.add_argument("--tss",              required=True,  help="TSS BED file (>=4 cols: chr start end gene)")
    ap.add_argument("--out",              required=True,  help="Output BEDPE file path")
    ap.add_argument("--score-col",        default="importance_x_abs_rho", help="Score column name in eRegulon file")
    ap.add_argument("--top-n",            type=int, default=300,           help="Keep top N links by score")
    ap.add_argument("--gene",             default=None,                    help="Filter to this target gene")
    ap.add_argument("--tf",               default=None,                    help="Filter to this TF")
    ap.add_argument("--region",           default=None,                    help="Filter to this locus (e.g. chr6:40963102-41363102)")
    ap.add_argument("--out-enhancers-bed",default=None,                    help="Optional: output BED of enhancer regions")
    ap.add_argument("--merge-enhancers",  action="store_true",             help="Deduplicate enhancer intervals before writing BED")
    ap.add_argument("--cmap",             default="Blues",                 help="Matplotlib colormap for arc coloring")
    ap.add_argument("--vmin",             type=float, default=None,        help="Min score for color scale")
    ap.add_argument("--vmax",             type=float, default=None,        help="Max score for color scale")
    ap.add_argument("--log-score",        action="store_true",             help="Apply log1p to scores before coloring/sorting")
    args = ap.parse_args()

    # Load inputs
    tss_map = load_tss(args.tss)
    df      = load_eregulon(args.ereg, tf=args.tf, gene=args.gene)

    # Parse enhancer regions
    parsed = df["Region"].apply(parse_region)
    df     = df.loc[parsed.notna()].copy()
    df[["chr1", "start1", "end1"]] = pd.DataFrame(parsed.dropna().tolist(), index=df.index)
    print(f"Loaded {len(df)} links after region parsing")

    # Map genes to TSS coordinates
    df = df.join(tss_map, on="Gene", how="inner")

    # Filter by locus
    if args.region:
        df = filter_by_locus(df, args.region)
        print(f"After locus filter: {len(df)} links")

    # Score and rank
    if args.score_col not in df.columns:
        raise ValueError(f"Score column '{args.score_col}' not found in eRegulon file")

    df["score"] = pd.to_numeric(df[args.score_col], errors="coerce").fillna(0.0)
    df["score_for_rank"] = np.log1p(df["score"]) if args.log_score else df["score"]
    df = df.sort_values("score_for_rank", ascending=False).head(args.top_n)

    # Write enhancer BED
    if args.out_enhancers_bed:
        enh = df[["chr1", "start1", "end1"]].copy()
        enh.columns = ["chr", "start", "end"]
        if args.merge_enhancers:
            enh = enh.drop_duplicates()
        enh["name"] = "enh"
        enh.to_csv(args.out_enhancers_bed, sep="\t", header=False, index=False)
        print(f"Wrote {len(enh)} enhancers -> {args.out_enhancers_bed}")

    # Map scores to RGB colors
    df["rgb"], used_vmin, used_vmax = score_to_rgb(
        df["score_for_rank"], args.cmap, args.vmin, args.vmax
    )

    # Write BEDPE
    out = df[["chr1", "start1", "end1", "chr2", "start2", "end2", "score", "rgb"]]
    out.to_csv(args.out, sep="\t", header=False, index=False)

    print(f"Wrote {len(out)} links -> {args.out}")
    print(f"Color scale ({'log1p' if args.log_score else 'linear'}): "
          f"vmin={used_vmin:.4g}, vmax={used_vmax:.4g}, cmap={args.cmap}")


if __name__ == "__main__":
    main()