from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

from .genesets import load_yaml, default_yaml
from .io import load_expr_matrix, read_table
from .core import score_modules, score_axes, summarise_used_genes
from .spatial import infer_spatial_keys, assign_rings, ring_summary

def _ensure_out(out: str) -> Path:
    p = Path(out)
    p.mkdir(parents=True, exist_ok=True)
    return p

def cmd_bulk(args) -> int:
    out = _ensure_out(args.out)
    gs = load_yaml(args.geneset or default_yaml("rna"))

    expr, _ = load_expr_matrix(args.expr, gene_col=args.gene_col, orientation=args.orientation)
    scores, used = score_modules(expr, gs.modules, min_genes=args.min_genes, zscore=not args.no_zscore)
    axes = score_axes(scores, gs.axes)

    out_df = pd.concat([scores, axes], axis=1)
    out_df.index.name = "sample_id"
    out_df.to_csv(out/"barriernet_bulk_scores.tsv", sep="\t")

    summarise_used_genes(used).to_csv(out/"barriernet_bulk_used_genes.tsv", sep="\t", index=False)
    panel_genes = sorted({g for gl in used.values() for g in gl})
    if panel_genes:
        expr.loc[:, [g for g in panel_genes if g in expr.columns]].to_csv(out/"barriernet_bulk_panel_matrix.tsv", sep="\t")

    print(f"[barriernet bulk] wrote outputs to: {out}")
    return 0

def cmd_ihc(args) -> int:
    out = _ensure_out(args.out)
    gs = load_yaml(args.geneset or default_yaml("ihc"))

    df = read_table(args.table)
    if args.sample_col not in df.columns:
        raise SystemExit(f"sample_col '{args.sample_col}' not found in columns")
    sample_col = args.sample_col

    drop_cols = set([sample_col] + (args.drop_cols or []))
    marker_cols = [c for c in df.columns if c not in drop_cols]

    X = df[marker_cols].copy()
    X.columns = [str(c).upper() for c in X.columns]
    for c in X.columns:
        X[c] = pd.to_numeric(X[c], errors="coerce")

    if args.agg == "median":
        Xg = pd.concat([df[[sample_col]], X], axis=1).groupby(sample_col).median(numeric_only=True)
    else:
        Xg = pd.concat([df[[sample_col]], X], axis=1).groupby(sample_col).mean(numeric_only=True)

    Xg.index = Xg.index.astype(str)
    scores, used = score_modules(Xg, gs.modules, min_genes=args.min_genes, zscore=not args.no_zscore)
    axes = score_axes(scores, gs.axes)

    out_df = pd.concat([scores, axes], axis=1)
    out_df.index.name = "sample_id"
    out_df.to_csv(out/"barriernet_ihc_scores.tsv", sep="\t")
    summarise_used_genes(used).to_csv(out/"barriernet_ihc_used_markers.tsv", sep="\t", index=False)

    print(f"[barriernet ihc] wrote outputs to: {out}")
    return 0

def _read_h5ad_any(path: str):
    try:
        import anndata as ad
        return ad.read_h5ad(path)
    except Exception:
        import scanpy as sc
        return sc.read_h5ad(path)

def cmd_spatial(args) -> int:
    out = _ensure_out(args.out)
    paths = []
    if args.h5ad:
        paths = [args.h5ad]
    else:
        import glob
        paths = sorted(glob.glob(args.h5ad_glob))

    if not paths:
        raise SystemExit("No h5ad files found.")

    for p in paths:
        adata = _read_h5ad_any(p)
        if "spatial" not in adata.obsm:
            raise SystemExit(f"{p}: missing obsm['spatial'].")

        xy = np.asarray(adata.obsm["spatial"]).astype(float)
        keys = infer_spatial_keys(adata.obs)

        ring_df = assign_rings(adata.obs, xy, keys.tls, peak_q=args.peak_q, peak_k=args.peak_k, far_min_dist=args.far_min_dist)

        stem = Path(p).name.replace(".h5ad","")
        spot = pd.concat([ring_df], axis=1)
        # Keep a minimal per-spot export that matches your Figure3 expectations
        spot.to_csv(out/f"{stem}_spot_rings.tsv", sep="\t")

        rs = ring_summary(adata.obs, ring_df, keys)
        rs.to_csv(out/f"{stem}_ring_enrichment.tsv", sep="\t", index=False)

        # Write back to h5ad
        for c in ring_df.columns:
            adata.obs[c] = ring_df[c].values
        out_h5 = out/f"{stem}.with_barriernet.h5ad"
        try:
            adata.write(out_h5)
        except Exception:
            import scanpy as sc
            sc.write(out_h5, adata)

        print(f"[barriernet spatial] processed: {stem}")

    print(f"[barriernet spatial] wrote outputs to: {out}")
    return 0

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="barriernet",
        description="Barrier/NetImmune portable scoring (bulk, IHC/IF) + Visium spatial ring/peak summaries."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    p_bulk = sub.add_parser("bulk", help="Compute Barrier/NetImmune axes from bulk RNA/proteomics expression matrix.")
    p_bulk.add_argument("--expr", required=True, help="Expression matrix TSV/CSV (samples×genes or genes×samples).")
    p_bulk.add_argument("--gene-col", default=None, help="Optional column name containing gene symbols.")
    p_bulk.add_argument("--orientation", default="auto", choices=["auto","samples_by_genes","genes_by_samples"])
    p_bulk.add_argument("--geneset", default=None, help="Path to geneset YAML (default: bundled RNA geneset).")
    p_bulk.add_argument("--min-genes", type=int, default=3)
    p_bulk.add_argument("--no-zscore", action="store_true")
    p_bulk.add_argument("--out", required=True)
    p_bulk.set_defaults(func=cmd_bulk)

    p_ihc = sub.add_parser("ihc", help="Compute Barrier/NetImmune axes from IHC/IF/protein panel table.")
    p_ihc.add_argument("--table", required=True, help="Long table with marker columns + sample_id column (ROI/cell rows allowed).")
    p_ihc.add_argument("--sample-col", required=True, help="Column name for sample/slide/patient ID.")
    p_ihc.add_argument("--drop-cols", nargs="*", default=None, help="Optional metadata columns to drop (besides sample_col).")
    p_ihc.add_argument("--agg", default="mean", choices=["mean","median"], help="Aggregate ROI/cell values to sample.")
    p_ihc.add_argument("--geneset", default=None, help="Path to geneset YAML (default: bundled IHC geneset).")
    p_ihc.add_argument("--min-genes", type=int, default=2)
    p_ihc.add_argument("--no-zscore", action="store_true")
    p_ihc.add_argument("--out", required=True)
    p_ihc.set_defaults(func=cmd_ihc)

    p_sp = sub.add_parser("spatial", help="Visium spatial: assign TLS-centered rings using hex-grid adjacency.")
    g = p_sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--h5ad", help="Single h5ad path.")
    g.add_argument("--h5ad-glob", help="Glob for multiple h5ad files.")
    p_sp.add_argument("--peak-q", type=float, default=0.95)
    p_sp.add_argument("--peak-k", type=int, default=9)
    p_sp.add_argument("--far-min-dist", type=int, default=4)
    p_sp.add_argument("--out", required=True)
    p_sp.set_defaults(func=cmd_spatial)

    # Backward compatible alias of spatial
    p_run = sub.add_parser("run", help="Alias of 'spatial' for backward compatibility.")
    g2 = p_run.add_mutually_exclusive_group(required=True)
    g2.add_argument("--h5ad", help="Single h5ad path.")
    g2.add_argument("--h5ad-glob", help="Glob for multiple h5ad files.")
    p_run.add_argument("--peak-q", type=float, default=0.95)
    p_run.add_argument("--peak-k", type=int, default=9)
    p_run.add_argument("--far-min-dist", type=int, default=4)
    p_run.add_argument("--out", required=True)
    p_run.set_defaults(func=cmd_spatial)

    return p

def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
