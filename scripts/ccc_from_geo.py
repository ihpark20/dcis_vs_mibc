"""Ligand-receptor communication, read off spatial proximity.

Two cells can only signal through a short-range ligand if they are close enough to touch
the same space. For each ligand-receptor pair on the panel, this asks whether cells
carrying the receptor sit next to cells carrying the ligand more often than the tissue's
own composition would give by chance:

    obs = ligand-positive neighbours of receptor-positive cells / all their neighbours
    exp = ligand-positive cells / all cells in the core
    enrichment = log2((obs + eps) / (exp + eps))

Dividing sums rather than averaging per-cell fractions keeps sparsely surrounded cells from
dominating, and `exp` is the core's own density, so a core full of ligand-positive cells is
not credited for proximity that comes for free. Positive enrichment means the two
populations are found together more than the core's composition predicts; it is evidence of
opportunity to signal, not of signalling.

Only pairs with both partners on the 374-gene panel can be tested — ten of them. A core
contributes a pair only if it holds at least 15 ligand-positive and 15 receptor-positive
cells; expression is a raw count above zero.

Cores are the replicates. Each pair is tested against zero enrichment with a one-sample
Wilcoxon corrected across pairs, and DCIS is compared with mDCIS by Mann-Whitney.

Usage:
    python scripts/use_published_subclusters.py
    python scripts/ccc_from_geo.py

Output:
    03.data_processed/ccc/lr_enrich.csv      enrichment per core and pair
    03.data_processed/ccc/lr_stats.csv       per pair, with the DCIS/mDCIS comparison
    03.data_processed/ccc/lr_celltypes.csv   which cell types carry each partner
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from analysis_exclusions import EXCLUDED_CORES
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent
CLUSTERED = ROOT / "03.data_processed/integrated_clusters.h5ad"
MAJOR = ROOT / "03.data_processed/subclustered/major_celltypes.csv"
OUT_DIR = ROOT / "03.data_processed/ccc"


RADIUS = 30.0
MIN_EXPR = 15
MIN_CORES = 5
MIN_CELLS_PER_CORE = 50
EPS = 1e-4

PAIRS = [
    ("CSF1", "CSF1R", "macrophage recruitment and survival"),
    ("CXCL13", "CXCR5", "B cell follicle, tertiary lymphoid structure"),
    ("CCL21", "CCR7", "T and dendritic cell homing"),
    ("CCL19", "CCR7", "T and dendritic cell homing"),
    ("CCL5", "CCR5", "T and NK recruitment"),
    ("CXCL12", "CXCR4", "stromal homing"),
    ("CD274", "PDCD1", "PD-L1 / PD-1 checkpoint"),
    ("CD80", "CTLA4", "CD80 / CTLA-4 checkpoint"),
    ("CD86", "CTLA4", "CD86 / CTLA-4 checkpoint"),
    ("AREG", "EGFR", "EGFR-driven tumor growth"),
]


def core_enrichment(pos, expressed, pairs):
    """Spatial enrichment of every testable pair within one core."""
    tree = cKDTree(pos)
    neighbours = tree.query_ball_point(pos, RADIUS)
    n_neighbours = np.array([len(n) - 1 for n in neighbours])

    out = []
    for ligand, receptor, _ in pairs:
        lig, rec = expressed[ligand], expressed[receptor]
        n_lig, n_rec = int(lig.sum()), int(rec.sum())
        if n_lig < MIN_EXPR or n_rec < MIN_EXPR:
            continue
        rec_idx = np.where(rec)[0]
        lig_neighbours = np.array(
            [lig[neighbours[i]].sum() - (1 if lig[i] else 0) for i in rec_idx]
        )
        totals = n_neighbours[rec_idx]
        usable = totals > 0
        if usable.sum() < MIN_EXPR:
            continue
        obs = lig_neighbours[usable].sum() / totals[usable].sum()
        exp = n_lig / len(pos)
        out.append(
            {
                "pair": f"{ligand}->{receptor}",
                "ligand": ligand,
                "receptor": receptor,
                "obs_frac": round(obs, 4),
                "exp_frac": round(exp, 4),
                "enrich": round(float(np.log2((obs + EPS) / (exp + EPS))), 3),
                "n_ligand": n_lig,
                "n_receptor": n_rec,
            }
        )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clustered", default=CLUSTERED, type=Path)
    ap.add_argument("--major", default=MAJOR, type=Path)
    ap.add_argument("--out-dir", default=OUT_DIR, type=Path)
    args = ap.parse_args()

    a = sc.read_h5ad(args.clustered)
    obs = a.obs.copy()
    obs["cell_id"] = a.obs_names
    obs["core"] = obs["slide"].astype(str) + "_" + obs["core_id"].astype(str)
    keep = (~obs["core"].isin(EXCLUDED_CORES)).to_numpy()

    counts = a.layers["counts"] if "counts" in a.layers else a.X
    counts = sp.csc_matrix(counts) if not sp.issparse(counts) else counts.tocsc()
    gene_index = {g: i for i, g in enumerate(a.var_names)}
    used = sorted({g for p in PAIRS for g in p[:2]})
    missing = [g for g in used if g not in gene_index]
    if missing:
        print(f"panel is missing {missing}; those pairs are skipped")
    expressed = {
        g: np.asarray(counts[np.where(keep)[0], gene_index[g]].todense()).ravel() > 0
        for g in used
        if g in gene_index
    }
    pairs = [p for p in PAIRS if p[0] in expressed and p[1] in expressed]

    obs = obs[keep]
    pos_all = obs[["x_centroid", "y_centroid"]].to_numpy()
    print(f"{len(obs):,} cells in {obs['core'].nunique()} cores, {len(pairs)} testable pairs")

    rows = []
    for core, idx in obs.groupby("core", observed=True).indices.items():
        if len(idx) < MIN_CELLS_PER_CORE:
            continue
        per_core = {g: v[idx] for g, v in expressed.items()}
        for rec in core_enrichment(pos_all[idx], per_core, pairs):
            rec["core"] = core
            rec["sample_group"] = obs["sample_group"].to_numpy()[idx][0]
            rows.append(rec)

    enrich = pd.DataFrame(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    enrich.to_csv(args.out_dir / "lr_enrich.csv", index=False)

    major = pd.read_csv(args.major, usecols=["cell_id", "major_celltype"])
    types = obs[["cell_id"]].merge(major, on="cell_id", how="left")["major_celltype"].to_numpy()
    ct_rows = []
    for ligand, receptor, _ in pairs:
        for side, gene in (("ligand", ligand), ("receptor", receptor)):
            share = pd.Series(types[expressed[gene]]).value_counts(normalize=True) * 100
            for celltype, pct in share.head(6).items():
                ct_rows.append(
                    {
                        "pair": f"{ligand}->{receptor}",
                        "side": side,
                        "gene": gene,
                        "major_celltype": celltype,
                        "pct": round(pct, 1),
                    }
                )
    celltypes = pd.DataFrame(ct_rows)
    celltypes.to_csv(args.out_dir / "lr_celltypes.csv", index=False)

    stats = []
    for ligand, receptor, description in pairs:
        pair = f"{ligand}->{receptor}"
        e = enrich[enrich["pair"] == pair]
        if len(e) < MIN_CORES:
            continue
        dcis = e.loc[e["sample_group"] == "dcis", "enrich"]
        mibc = e.loc[e["sample_group"] == "mibc", "enrich"]
        stats.append(
            {
                "pair": pair,
                "description": description,
                "n_cores": len(e),
                "median_enrich": round(e["enrich"].median(), 3),
                "dcis_median": round(dcis.median(), 3) if len(dcis) else np.nan,
                "mibc_median": round(mibc.median(), 3) if len(mibc) else np.nan,
                "p_enrich": wilcoxon(e["enrich"]).pvalue if (e["enrich"] != 0).any() else 1.0,
                "p_dcis_vs_mibc": (
                    mannwhitneyu(dcis, mibc).pvalue if min(len(dcis), len(mibc)) >= 3 else np.nan
                ),
            }
        )
    stats = pd.DataFrame(stats)
    stats["fdr_enrich"] = multipletests(stats["p_enrich"], method="fdr_bh")[1]
    stats = stats.sort_values("median_enrich", ascending=False)
    stats.to_csv(args.out_dir / "lr_stats.csv", index=False)

    print("\n" + stats.round(4).to_string(index=False))
    print(f"\nwrote {args.out_dir}")


if __name__ == "__main__":
    main()
