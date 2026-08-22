"""CL5 + CL8 통합 B/Plasma 서브클러스터링.

CL5 (B cell 우세, 52,097) + CL8 (plasma cell 우세, 16,426)을 합쳐 재클러스터링.
B/plasma 마커 하나라도 양성인 세포만 필터 후, B/plasma focused feature로
PCA→Harmony→Leiden. 각 서브클러스터 subtype/마커/프로그램 + 코어 dominance.

전략:
  1. CL5+CL8 subset 결합
  2. B/plasma 마커(>0) 하나라도 양성 필터
  3. B/plasma focused feature로 PCA→Harmony(slide)→Leiden
  4. 프로그램 스코어 + parent(CL5/CL8) 분포 + 코어 dominance(Top1/2/3)

출력:
    03.data_processed/cl5_8_bplasma_subclustered.h5ad
    04.tables/cl5_8_bplasma_marker_genes.csv
    04.tables/cl5_8_bplasma_program_scores.csv
    04.tables/cl5_8_bplasma_summary.csv          (subtype/마커/프로그램/dominance)
    04.tables/cl5_8_bplasma_parent_distribution.csv
    05.figures/cl5_8_bplasma_umap_subclusters.png
    05.figures/cl5_8_bplasma_umap_program_scores.png
    05.figures/cl5_8_bplasma_dotplot_markers.png
    05.figures/cl5_8_bplasma_program_score_heatmap.png
"""

import sys

sys.modules["tensorflow"] = None  # noqa: E402

from pathlib import Path  # noqa: E402

import anndata  # noqa: E402
import harmonypy as hm  # noqa: E402
import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402

matplotlib.use("Agg")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

# Published verbatim except for the four path constants below, which point at this
# repository instead of the original project, and at the object prepare_pool_input.py
# builds. Everything that decides a result — filters, features, parameters, seeds, the
# programme scoring and the subtype rule — is the code that produced the paper.

ROOT = Path(__file__).resolve().parent.parent
PROJ = ROOT
INPUT_H5AD = PROJ / "03.data_processed/pool_input.h5ad"
OUT_H5AD = PROJ / "03.data_processed/subclustered/cl5_8_bplasma_subclustered.h5ad"
FIG_DIR = PROJ / "03.data_processed/subclustered/figures"
TBL_DIR = PROJ / "03.data_processed/subclustered"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 1

TARGET_CLUSTERS = ["5", "8"]
RESOLUTION = 0.4
RANDOM_STATE = 42
N_PCS = 15
LABEL = "cl5_8_bplasma"
QC_FAIL = {("TMA1", "D4"), ("TMA1", "D8")}

# B/plasma 마커 하나라도 양성 → 유지
FILTER_GENES = [
    "MS4A1",
    "CD79A",
    "CD79B",
    "CD19",
    "BANK1",
    "TCL1A",
    "CD27",
    "MZB1",
    "TNFRSF17",
    "CD38",
    "PRDM1",
    "IRF4",
    "IGHG1",
    "IGHM",
]

FEATURE_GENES = [
    # pan-B
    "MS4A1",
    "CD79A",
    "CD79B",
    "CD19",
    "BANK1",
    # naive/transitional
    "TCL1A",
    "IGHD",
    "FCER2",
    "SELL",
    "CR2",
    # memory / activation
    "CD27",
    "CD83",
    "CD40",
    # germinal center
    "BCL6",
    "CXCR5",
    "AICDA",
    "RGS13",
    # atypical / age-associated B
    "FCRL4",
    "ITGAX",
    "TBX21",
    # plasma / plasmablast
    "MZB1",
    "TNFRSF17",
    "PRDM1",
    "IRF4",
    "CD38",
    "XBP1",
    "SLAMF7",
    # Ig isotype
    "IGHG1",
    "IGHM",
    "IGHA1",
    "IGKC",
    "JCHAIN",
    # proliferation
    "MKI67",
    "TOP2A",
    "CENPF",
    "HMGA1",
    # antigen presentation
    "HLA-DRA",
    "CD74",
    "CXCR4",
]

PROGRAM_MARKERS = {
    "pan_B": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "naive_B": ["TCL1A", "IGHD", "FCER2", "SELL", "CR2"],
    "memory_B": ["CD27", "CD83"],
    "germinal_center": ["BCL6", "CXCR5", "AICDA", "RGS13"],
    "atypical_B": ["FCRL4", "ITGAX", "TBX21"],
    "plasma": ["MZB1", "TNFRSF17", "PRDM1", "IRF4", "CD38", "XBP1", "SLAMF7"],
    "igG": ["IGHG1"],
    "igM": ["IGHM"],
    "proliferating": ["MKI67", "TOP2A", "CENPF", "HMGA1"],
    "antigen_presentation": ["HLA-DRA", "CD74"],
}

# subtype 라벨 결정에 쓸 프로그램 (pan_B/isotype/AP 제외)
SUBTYPE_PROGS = ["naive_B", "memory_B", "germinal_center", "atypical_B", "plasma", "proliferating"]


def main():
    print("=" * 60)
    print(f"1. CL{'/'.join(TARGET_CLUSTERS)} subset 결합 + filter")
    print("=" * 60)
    adata_full = anndata.read_h5ad(INPUT_H5AD)
    sub = adata_full[adata_full.obs["leiden"].isin(TARGET_CLUSTERS)].copy()
    # QC FAIL core 제거
    for s, c in QC_FAIL:
        sub = sub[~((sub.obs["slide"] == s) & (sub.obs["core_id"] == c))].copy()
    print("  결합 전:")
    print(sub.obs["leiden"].value_counts().sort_index().to_string())
    n_before = sub.n_obs

    counts = sub.layers["counts"]
    filter_idx = [sub.var_names.get_loc(g) for g in FILTER_GENES if g in sub.var_names]
    print(f"  filter genes({len(filter_idx)}): {[g for g in FILTER_GENES if g in sub.var_names]}")
    sub_counts = counts[:, filter_idx].toarray() if sp.issparse(counts) else counts[:, filter_idx]
    is_pass = (sub_counts > 0).any(axis=1)
    sub = sub[is_pass].copy()
    print(f"  filter > 0: {sub.n_obs:,} / {n_before:,} ({sub.n_obs / n_before * 100:.1f}%)")
    print(f"  병리: {sub.obs['pathology'].value_counts().to_dict()}")

    sub.obs["leiden_parent"] = sub.obs["leiden"].astype(str).map(lambda x: f"CL{x}")
    sub.obs = sub.obs.drop(columns=["leiden"])

    print("\n2. 정규화 + log1p")
    sub.X = sub.layers["counts"].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sub.layers["lognorm"] = sub.X.copy()

    print("\n3. Feature 선택")
    feature_genes = [g for g in FEATURE_GENES if g in sub.var_names]
    print(f"  사용 마커: {len(feature_genes)}/{len(FEATURE_GENES)}")
    sub.var["highly_variable"] = sub.var_names.isin(feature_genes)

    print("\n4. Scale + PCA")
    sc.pp.scale(sub, max_value=10)
    n_comps = min(N_PCS, len(feature_genes) - 1)
    sc.tl.pca(sub, n_comps=n_comps, use_highly_variable=True)

    print("\n5. Harmony (slide)")
    harmony_out = hm.run_harmony(
        sub.obsm["X_pca"].copy(),
        sub.obs[["slide"]].copy(),
        vars_use=["slide"],
        random_state=RANDOM_STATE,
        max_iter_harmony=30,
    )
    sub.obsm["X_pca_harmony"] = np.array(harmony_out.Z_corr)

    print(f"\n6. Neighbors + UMAP + Leiden (res={RESOLUTION})")
    sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=n_comps)
    sc.tl.umap(sub, random_state=RANDOM_STATE)
    sc.tl.leiden(sub, resolution=RESOLUTION, random_state=RANDOM_STATE, key_added="leiden_sub")
    n_sub = sub.obs["leiden_sub"].nunique()
    print(f"  → 서브클러스터 수: {n_sub}")
    print(sub.obs["leiden_sub"].value_counts().sort_index().to_string())

    print("\n7. 프로그램 스코어")
    sub_score = sub.copy()
    sub_score.X = sub_score.layers["lognorm"].copy()
    program_keys = []
    for prog, markers in PROGRAM_MARKERS.items():
        present = [g for g in markers if g in sub_score.var_names]
        if not present:
            continue
        key = f"score_{prog}"
        sc.tl.score_genes(sub_score, gene_list=present, score_name=key, random_state=RANDOM_STATE)
        sub.obs[key] = sub_score.obs[key].values
        program_keys.append(key)
        print(f"  {prog}: {len(present)} markers")

    print("\n8. 서브클러스터별 마커 (Wilcoxon)")
    sc.tl.rank_genes_groups(
        sub_score,
        groupby="leiden_sub",
        method="wilcoxon",
        key_added="rank_sub",
        n_genes=50,
        use_raw=False,
    )
    sub_cats = sorted(sub.obs["leiden_sub"].unique(), key=int)
    rec = []
    for cl in sub_cats:
        names = sub_score.uns["rank_sub"]["names"][cl]
        scores = sub_score.uns["rank_sub"]["scores"][cl]
        pvals = sub_score.uns["rank_sub"]["pvals_adj"][cl]
        lfcs = sub_score.uns["rank_sub"]["logfoldchanges"][cl]
        for r, (g, s, p, lfc) in enumerate(zip(names, scores, pvals, lfcs, strict=False), start=1):
            rec.append(
                {
                    "leiden_sub": cl,
                    "rank": r,
                    "gene": g,
                    "score": round(float(s), 3),
                    "logfoldchange": round(float(lfc), 3),
                    "pval_adj": float(p),
                }
            )
    marker_df = pd.DataFrame(rec)
    marker_df.to_csv(TBL_DIR / f"{LABEL}_marker_genes.csv", index=False)

    print("\n9. 프로그램 평균")
    score_mat = sub.obs.groupby("leiden_sub", observed=True)[program_keys].mean().loc[sub_cats]
    score_mat.round(4).to_csv(TBL_DIR / f"{LABEL}_program_scores.csv")
    print(score_mat.round(2).to_string())

    print("\n10. parent(CL5/CL8) 분포")
    parent_dist = (
        sub.obs.groupby(["leiden_sub", "leiden_parent"], observed=True)
        .size()
        .unstack(fill_value=0)
        .loc[sub_cats]
    )
    parent_pct = parent_dist.div(parent_dist.sum(axis=1), axis=0) * 100
    parent_combined = pd.concat(
        [parent_dist.add_prefix("n_"), parent_pct.add_prefix("pct_").round(1)], axis=1
    )
    parent_combined.to_csv(TBL_DIR / f"{LABEL}_parent_distribution.csv")

    print("\n11. 서브클러스터 요약 + 코어 dominance(Top1/2/3)")
    sub.obs["core"] = sub.obs["slide"].astype(str) + "_" + sub.obs["core_id"].astype(str)
    core_path = sub.obs.drop_duplicates("core").set_index("core")["pathology"].to_dict()
    summary_rows = []
    for cl in sub_cats:
        cells = sub.obs[sub.obs["leiden_sub"] == cl]
        n = len(cells)
        path_pct = cells["pathology"].value_counts(normalize=True) * 100
        parent_pct_cl = cells["leiden_parent"].value_counts(normalize=True) * 100
        ct_dist = cells["celltype_final"].value_counts(normalize=True) * 100
        top_markers = marker_df[marker_df["leiden_sub"] == cl].head(8)["gene"].tolist()
        prog_means = score_mat.loc[cl].sort_values(ascending=False)
        top_prog = ", ".join(
            [f"{k.replace('score_', '')}({v:.2f})" for k, v in prog_means.head(3).items()]
        )
        subtype_means = {p: score_mat.loc[cl].get(f"score_{p}", -np.inf) for p in SUBTYPE_PROGS}
        subtype_label = max(subtype_means, key=subtype_means.get)
        # 코어 dominance
        cc = cells["core"].value_counts()
        share = (cc / cc.sum() * 100).values
        cum = np.cumsum(share)
        tops = cc.index.tolist()
        hhi = float(np.sum((share / 100) ** 2))
        summary_rows.append(
            {
                "leiden_sub": cl,
                "n_cells": n,
                "subtype": subtype_label,
                "pct_CL5": round(parent_pct_cl.get("CL5", 0), 1),
                "pct_CL8": round(parent_pct_cl.get("CL8", 0), 1),
                "pct_mibc": round(path_pct.get("mibc", 0), 1),
                "pct_dcis": round(path_pct.get("dcis", 0), 1),
                "dominant_prev_celltype": ct_dist.index[0] if len(ct_dist) else "unknown",
                "dominant_prev_pct": round(ct_dist.iloc[0], 1) if len(ct_dist) else 0.0,
                "top1_core": f"{tops[0]}({core_path.get(tops[0], '?')})",
                "top1_pct": round(share[0], 1),
                "top2_cum_pct": round(cum[1], 1) if len(cum) > 1 else round(cum[-1], 1),
                "top3_cum_pct": round(cum[2], 1) if len(cum) > 2 else round(cum[-1], 1),
                "eff_cores": round(1 / hhi, 1),
                "top3_cores": "; ".join(
                    f"{c}({core_path.get(c, '?')}, {share[i]:.0f}%)" for i, c in enumerate(tops[:3])
                ),
                "top_programs": top_prog,
                "top8_markers": ", ".join(top_markers),
            }
        )
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(TBL_DIR / f"{LABEL}_summary.csv", index=False)
    print(
        summary_df[
            [
                "leiden_sub",
                "n_cells",
                "subtype",
                "pct_CL5",
                "pct_CL8",
                "pct_mibc",
                "top1_pct",
                "top3_cum_pct",
                "eff_cores",
                "top8_markers",
            ]
        ].to_string(index=False)
    )

    print("\n12. Figures")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    sc.pl.umap(sub, color="leiden_sub", legend_loc="on data", ax=axes[0], show=False, frameon=False)
    axes[0].set_title(f"CL5+CL8 sub-Leiden (res={RESOLUTION}, n={n_sub})")
    sc.pl.umap(
        sub,
        color="leiden_parent",
        ax=axes[1],
        show=False,
        frameon=False,
        palette={"CL5": "#1f77b4", "CL8": "#d62728"},
    )
    axes[1].set_title("Parent (CL5=B / CL8=plasma)")
    sc.pl.umap(sub, color="pathology", ax=axes[2], show=False, frameon=False)
    axes[2].set_title("Pathology")
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{LABEL}_umap_subclusters.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    n_prog = len(program_keys)
    n_col = 4
    n_row = int(np.ceil(n_prog / n_col))
    fig, axes = plt.subplots(n_row, n_col, figsize=(5 * n_col, 4.5 * n_row))
    axes = np.array(axes).reshape(-1)
    for i, key in enumerate(program_keys):
        sc.pl.umap(
            sub,
            color=key,
            ax=axes[i],
            show=False,
            frameon=False,
            cmap="viridis",
            vmin="p2",
            vmax="p98",
        )
        axes[i].set_title(key.replace("score_", ""))
    for j in range(len(program_keys), len(axes)):
        axes[j].axis("off")
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{LABEL}_umap_program_scores.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig_dot = sc.pl.dotplot(
        sub_score,
        var_names=feature_genes,
        groupby="leiden_sub",
        standard_scale="var",
        show=False,
        return_fig=True,
        figsize=(max(12, len(feature_genes) * 0.4), max(4, n_sub * 0.45)),
    )
    fig_dot.savefig(FIG_DIR / f"{LABEL}_dotplot_markers.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    score_z = ((score_mat - score_mat.mean()) / score_mat.std()).clip(-3, 3)
    fig, ax = plt.subplots(figsize=(max(8, len(program_keys) * 0.7), max(4, n_sub * 0.4)))
    im = ax.imshow(score_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(program_keys)))
    ax.set_xticklabels([k.replace("score_", "") for k in program_keys], rotation=45, ha="right")
    ax.set_yticks(range(n_sub))
    ax.set_yticklabels([f"sub{c}" for c in sub_cats])
    plt.colorbar(im, ax=ax, label="z-score", shrink=0.7)
    for i in range(n_sub):
        for j in range(len(program_keys)):
            ax.text(j, i, f"{score_mat.values[i, j]:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("CL5+CL8 × 프로그램 스코어 (z-score)", fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{LABEL}_program_score_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("\n13. h5ad 저장")
    sub.write_h5ad(OUT_H5AD)
    print(f"  → {OUT_H5AD}\n완료!")


if __name__ == "__main__":
    main()
