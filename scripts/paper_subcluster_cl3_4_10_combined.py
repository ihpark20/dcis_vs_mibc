"""
CL3 + CL4 + CL10 통합 재클러스터링

세 클러스터 모두 tumor 관련:
  - CL3 (55,809): HER2+ tumor 풀
  - CL4 (52,806): tumor + myoepithelial 혼재
  - CL10 (8,274): HER2- MIBC tumor 풀

합쳐서 재클러스터링하면 HER2+ vs HER2-, tumor vs myoepithelial 분리 양상을
한 공간에서 비교 가능. 또 클러스터 간 공통 state(예: luminal A, EMT,
proliferating)가 합쳐지는지도 볼 수 있음.

전략:
  1. CL3+CL4+CL10 subset 결합
  2. Epithelial+myoep filter (KRT8/EPCAM/CDH1/KRT7/TACSTD2/KRT5/KRT14/ACTA2/MYH11)
  3. Tumor+myoep focused feature(~37)로 PCA→Harmony→Leiden
  4. 프로그램 스코어 + parent cluster 분포 추적

출력:
    03.data_processed/cl3_4_10_combined_subclustered.h5ad
    04.tables/cl3_4_10_combined_marker_genes.csv
    04.tables/cl3_4_10_combined_program_scores.csv
    04.tables/cl3_4_10_combined_summary.csv
    04.tables/cl3_4_10_combined_parent_distribution.csv
    05.figures/cl3_4_10_combined_umap_subclusters.png
    05.figures/cl3_4_10_combined_umap_program_scores.png
    05.figures/cl3_4_10_combined_dotplot_markers.png
    05.figures/cl3_4_10_combined_program_score_heatmap.png
    05.figures/cl3_4_10_combined_umap_parent.png
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
OUT_H5AD = PROJ / "03.data_processed/subclustered/cl3_4_10_combined_subclustered.h5ad"
FIG_DIR = PROJ / "03.data_processed/subclustered/figures"
TBL_DIR = PROJ / "03.data_processed/subclustered"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 2

TARGET_CLUSTERS = ["3", "4", "10"]
RESOLUTION = 0.3
RANDOM_STATE = 42
N_PCS = 15
LABEL = "cl3_4_10_combined"

# Epithelial + myoepithelial filter (CL4가 myoep 포함하므로 둘 다)
FILTER_GENES = ["KRT8", "EPCAM", "CDH1", "KRT7", "TACSTD2", "KRT5", "KRT14", "ACTA2", "MYH11"]

FEATURE_GENES = [
    # pan-epithelial
    "KRT8",
    "EPCAM",
    "CDH1",
    "KRT7",
    "TACSTD2",
    # luminal mature
    "ESR1",
    "PGR",
    "FOXA1",
    "GATA3",
    "ANKRD30A",
    "AGR3",
    "MLPH",
    "BCL2",
    "AREG",
    "CCND1",
    # luminal progenitor
    "KIT",
    "ALDH1A3",
    "KRT23",
    "SFRP1",
    # HER2
    "ERBB2",
    # basal / triple-negative
    "KRT5",
    "KRT14",
    "KRT6B",
    "EGFR",
    # myoepithelial-specific
    "ACTA2",
    "MYH11",
    "MYLK",
    "ACTG2",
    "OXTR",
    # proliferation
    "MKI67",
    "TOP2A",
    "CENPF",
    "HMGA1",
    # EMT / invasion
    "ZEB1",
    "ZEB2",
    "SNAI1",
    "FN1",
    "S100A4",
    "POSTN",
    # cancer-specific differentiation
    "CEACAM6",
    "TPD52",
    "ERBB3",
    # stress / cell cycle inhibitor
    "CDKN1A",
]

PROGRAM_MARKERS = {
    "luminal_mature": [
        "ESR1",
        "PGR",
        "FOXA1",
        "GATA3",
        "ANKRD30A",
        "AGR3",
        "MLPH",
        "BCL2",
        "AREG",
    ],
    "luminal_progenitor": ["KIT", "ALDH1A3", "KRT23", "SFRP1"],
    "her2": ["ERBB2"],
    "basal_triple_neg": ["KRT5", "KRT14", "KRT6B", "EGFR"],
    "myoepithelial": ["ACTA2", "MYH11", "MYLK", "ACTG2", "OXTR"],
    "pan_epithelial": ["KRT8", "EPCAM", "CDH1", "KRT7", "TACSTD2"],
    "proliferating": ["MKI67", "TOP2A", "CENPF", "HMGA1"],
    "emt": ["ZEB1", "ZEB2", "SNAI1", "FN1", "S100A4", "POSTN"],
    "cancer_specific": ["CEACAM6", "TPD52", "AGR3", "MLPH"],
    "stress_g1arrest": ["CDKN1A"],
}


def main():
    print("=" * 60)
    print(f"1. CL{'/'.join(TARGET_CLUSTERS)} subset 결합 + filter")
    print("=" * 60)
    adata_full = anndata.read_h5ad(INPUT_H5AD)
    sub = adata_full[adata_full.obs["leiden"].isin(TARGET_CLUSTERS)].copy()
    print("  결합 전 원본:")
    print(sub.obs["leiden"].value_counts().sort_index().to_string())
    n_before = sub.n_obs
    print(f"  total: {n_before:,} cells")

    counts = sub.layers["counts"]
    filter_idx = [sub.var_names.get_loc(g) for g in FILTER_GENES if g in sub.var_names]
    print(f"  filter genes: {[g for g in FILTER_GENES if g in sub.var_names]}")
    sub_counts = counts[:, filter_idx].toarray() if sp.issparse(counts) else counts[:, filter_idx]
    is_pass = (sub_counts > 0).any(axis=1)
    sub = sub[is_pass].copy()
    print(f"  filter > 0: {sub.n_obs:,} / {n_before:,} ({sub.n_obs / n_before * 100:.1f}%)")
    print(f"  슬라이드: {sub.obs['slide'].value_counts().to_dict()}")
    print(f"  병리: {sub.obs['pathology'].value_counts().to_dict()}")

    # parent 클러스터 라벨 보존
    sub.obs["leiden_parent"] = sub.obs["leiden"].astype(str).map(lambda x: f"CL{x}")
    sub.obs = sub.obs.drop(columns=["leiden"])

    # ── 정규화 ─────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("2. 정규화 + log1p")
    print("=" * 60)
    sub.X = sub.layers["counts"].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sub.layers["lognorm"] = sub.X.copy()

    # ── feature 선택 ───────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("3. Feature 선택 (epi + myoep + tumor states)")
    print("=" * 60)
    feature_genes = [g for g in FEATURE_GENES if g in sub.var_names]
    print(f"  사용 마커: {len(feature_genes)}/{len(FEATURE_GENES)}")
    sub.var["highly_variable"] = sub.var_names.isin(feature_genes)

    # ── Scale + PCA ────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print(f"4. Scale + PCA (n_comps={N_PCS})")
    print("=" * 60)
    sc.pp.scale(sub, max_value=10)
    n_comps = min(N_PCS, len(feature_genes) - 1)
    sc.tl.pca(sub, n_comps=n_comps, use_highly_variable=True)

    # ── Harmony ────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("5. Harmony 배치 보정 (slide)")
    print("=" * 60)
    pca_mat = sub.obsm["X_pca"].copy()
    harmony_out = hm.run_harmony(
        pca_mat,
        sub.obs[["slide"]].copy(),
        vars_use=["slide"],
        random_state=RANDOM_STATE,
        max_iter_harmony=30,
    )
    sub.obsm["X_pca_harmony"] = np.array(harmony_out.Z_corr)

    # ── Neighbors + UMAP + Leiden ──────────────────────────────────────────
    print("\n" + "=" * 60)
    print(f"6. Neighbors + UMAP + Leiden (res={RESOLUTION})")
    print("=" * 60)
    sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=n_comps)
    sc.tl.umap(sub, random_state=RANDOM_STATE)
    sc.tl.leiden(sub, resolution=RESOLUTION, random_state=RANDOM_STATE, key_added="leiden_sub")
    n_sub = sub.obs["leiden_sub"].nunique()
    print(f"  → 서브클러스터 수: {n_sub}")
    print(sub.obs["leiden_sub"].value_counts().sort_index().to_string())

    # ── 프로그램 스코어 ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("7. 프로그램 스코어")
    print("=" * 60)
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
        print(f"  {prog}: {len(present)} markers → {key}")

    # ── 마커 ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("8. 서브클러스터별 마커 (Wilcoxon)")
    print("=" * 60)
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
    marker_out = TBL_DIR / f"{LABEL}_marker_genes.csv"
    marker_df.to_csv(marker_out, index=False)
    print(f"  → {marker_out}")

    # ── 프로그램 평균 ──────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("9. 서브클러스터별 프로그램 스코어 평균")
    print("=" * 60)
    score_mat = sub.obs.groupby("leiden_sub", observed=True)[program_keys].mean().loc[sub_cats]
    score_out = TBL_DIR / f"{LABEL}_program_scores.csv"
    score_mat.round(4).to_csv(score_out)
    print(score_mat.round(2).to_string())

    # ── 부모 클러스터 분포 ─────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("10. 서브클러스터 × parent (CL3/CL4/CL10) 분포")
    print("=" * 60)
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
    parent_out = TBL_DIR / f"{LABEL}_parent_distribution.csv"
    parent_combined.to_csv(parent_out)
    print(f"  → {parent_out}")
    print(parent_combined.to_string())

    # ── 요약 ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("11. 서브클러스터 요약")
    print("=" * 60)
    summary_rows = []
    for cl in sub_cats:
        cells = sub.obs[sub.obs["leiden_sub"] == cl]
        n = len(cells)
        slide_pct = cells["slide"].value_counts(normalize=True) * 100
        path_pct = cells["pathology"].value_counts(normalize=True) * 100
        parent_pct_cl = cells["leiden_parent"].value_counts(normalize=True) * 100
        ct_dist = cells["celltype_final"].value_counts(normalize=True) * 100
        top_markers = marker_df[marker_df["leiden_sub"] == cl].head(8)["gene"].tolist()
        prog_means = score_mat.loc[cl].sort_values(ascending=False)
        top_prog = ", ".join(
            [f"{k.replace('score_', '')}({v:.2f})" for k, v in prog_means.head(3).items()]
        )
        # 가장 우세 프로그램 기준 라벨
        subtype_progs = [
            "luminal_mature",
            "luminal_progenitor",
            "her2",
            "basal_triple_neg",
            "myoepithelial",
            "emt",
            "proliferating",
        ]
        subtype_means = {p: score_mat.loc[cl].get(f"score_{p}", -np.inf) for p in subtype_progs}
        subtype_label = max(subtype_means, key=subtype_means.get)
        summary_rows.append(
            {
                "leiden_sub": cl,
                "n_cells": n,
                "subtype": subtype_label,
                "pct_CL3": round(parent_pct_cl.get("CL3", 0), 1),
                "pct_CL4": round(parent_pct_cl.get("CL4", 0), 1),
                "pct_CL10": round(parent_pct_cl.get("CL10", 0), 1),
                "pct_TMA1": round(slide_pct.get("TMA1", 0), 1),
                "pct_TMA2": round(slide_pct.get("TMA2", 0), 1),
                "pct_mibc": round(path_pct.get("mibc", 0), 1),
                "pct_dcis": round(path_pct.get("dcis", 0), 1),
                "dominant_prev_celltype": ct_dist.index[0] if len(ct_dist) > 0 else "unknown",
                "dominant_prev_pct": round(ct_dist.iloc[0], 1) if len(ct_dist) > 0 else 0.0,
                "top_programs": top_prog,
                "top8_markers": ", ".join(top_markers),
            }
        )
    summary_df = pd.DataFrame(summary_rows)
    summary_out = TBL_DIR / f"{LABEL}_summary.csv"
    summary_df.to_csv(summary_out, index=False)
    print(summary_df.to_string(index=False))

    # ── Figures ────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("12. Figures")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    sc.pl.umap(sub, color="leiden_sub", legend_loc="on data", ax=axes[0], show=False, frameon=False)
    axes[0].set_title(f"CL3+CL4+CL10 sub-Leiden (res={RESOLUTION}, n={n_sub})")
    sc.pl.umap(sub, color="slide", ax=axes[1], show=False, frameon=False)
    axes[1].set_title("Slide")
    sc.pl.umap(sub, color="pathology", ax=axes[2], show=False, frameon=False)
    axes[2].set_title("Pathology")
    plt.tight_layout()
    out = FIG_DIR / f"{LABEL}_umap_subclusters.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")

    # parent UMAP — 부모 클러스터별 색칠
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(
        sub,
        color="leiden_parent",
        ax=ax,
        show=False,
        frameon=False,
        palette={"CL3": "#1f77b4", "CL4": "#2ca02c", "CL10": "#d62728"},
    )
    ax.set_title("CL3+CL4+CL10 통합 UMAP (parent 클러스터)")
    plt.tight_layout()
    out = FIG_DIR / f"{LABEL}_umap_parent.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")

    # program score UMAPs
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
    out = FIG_DIR / f"{LABEL}_umap_program_scores.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")

    fig_dot = sc.pl.dotplot(
        sub_score,
        var_names=feature_genes,
        groupby="leiden_sub",
        standard_scale="var",
        show=False,
        return_fig=True,
        title="CL3+CL4+CL10 통합 서브클러스터 × 마커",
        figsize=(max(12, len(feature_genes) * 0.4), max(4, n_sub * 0.45)),
    )
    out = FIG_DIR / f"{LABEL}_dotplot_markers.png"
    fig_dot.savefig(out, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  → {out}")

    score_z = (score_mat - score_mat.mean()) / score_mat.std()
    score_z = score_z.clip(-3, 3)
    fig, ax = plt.subplots(figsize=(max(8, len(program_keys) * 0.7), max(4, n_sub * 0.4)))
    im = ax.imshow(score_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(program_keys)))
    ax.set_xticklabels([k.replace("score_", "") for k in program_keys], rotation=45, ha="right")
    ax.set_yticks(range(n_sub))
    ax.set_yticklabels([f"sub{c}" for c in sub_cats])
    ax.set_title("CL3+CL4+CL10 통합 × 프로그램 스코어 (z-score)", fontweight="bold")
    plt.colorbar(im, ax=ax, label="z-score", shrink=0.7)
    for i in range(n_sub):
        for j in range(len(program_keys)):
            ax.text(
                j,
                i,
                f"{score_mat.values[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=7,
                color="black",
            )
    plt.tight_layout()
    out = FIG_DIR / f"{LABEL}_program_score_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")

    # ── 저장 ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("13. h5ad 저장")
    print("=" * 60)
    sub.write_h5ad(OUT_H5AD)
    print(f"  → {OUT_H5AD}")
    print("\n완료!")


if __name__ == "__main__":
    main()
