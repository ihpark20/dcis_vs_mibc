"""
CL1 (endothelial cell dominant) 서브클러스터링

CL1 (58,750 cells, dominant=endothelial cell, 93.8% stroma label)을
CL0/CL2/CL4와 동일 전략으로 분석:
  1. CL1 subset
  2. Endothelial filter (PECAM1/VWF/KDR any > 0)
  3. EC focused feature(~21)로 PCA→Harmony→Leiden
  4. 프로그램 스코어:
     - pan_endothelial, arterial, capillary, tip_sprouting,
       hev_like, activated_inflam, pericyte_mural, endmt, proliferating
  5. Tables, figures

Note: Panel에 venous/lymphatic 마커가 거의 없어서 arterial-venous 분류는 제한적.

출력:
    03.data_processed/cl1_endothelial_subclustered.h5ad
    04.tables/cl1_marker_genes.csv
    04.tables/cl1_program_scores.csv
    04.tables/cl1_summary.csv
    05.figures/cl1_umap_subclusters.png
    05.figures/cl1_umap_program_scores.png
    05.figures/cl1_dotplot_markers.png
    05.figures/cl1_program_score_heatmap.png
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
OUT_H5AD = PROJ / "03.data_processed/subclustered/cl1_endothelial_subclustered.h5ad"
FIG_DIR = PROJ / "03.data_processed/subclustered/figures"
TBL_DIR = PROJ / "03.data_processed/subclustered"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 2

TARGET_CLUSTER = "1"
RESOLUTION = 0.3
RANDOM_STATE = 42
N_PCS = 15

EC_FILTER_GENES = ["PECAM1", "VWF", "KDR"]

FEATURE_GENES = [
    # pan-endothelial
    "PECAM1",
    "VWF",
    "KDR",
    # arterial (제한적)
    "SOX17",
    # capillary
    "AQP1",
    "RAMP2",
    # tip cell / sprouting angiogenesis
    "ESM1",
    "ANGPT2",
    "CXCR4",
    # high endothelial venule (HEV-like)
    "CCL19",
    "CCL21",
    # activated/inflammatory
    "CXCL12",
    # pericyte / mural cell distinction
    "ACTA2",
    "MYH11",
    "PDGFRB",
    # EndMT
    "ZEB1",
    "ZEB2",
    "SNAI1",
    "FN1",
    "S100A4",
    # proliferation
    "MKI67",
    "TOP2A",
]

PROGRAM_MARKERS = {
    "pan_endothelial": ["PECAM1", "VWF", "KDR"],
    "arterial": ["SOX17"],
    "capillary": ["AQP1", "RAMP2"],
    "tip_sprouting": ["ESM1", "ANGPT2", "CXCR4"],
    "hev_like": ["CCL19", "CCL21"],
    "activated_inflam": ["CXCL12"],
    "pericyte_mural": ["ACTA2", "MYH11", "PDGFRB"],
    "endmt": ["ZEB1", "ZEB2", "SNAI1", "FN1", "S100A4"],
    "proliferating": ["MKI67", "TOP2A"],
}


def main():
    print("=" * 60)
    print(f"1. CL{TARGET_CLUSTER} subset 추출 + endothelial filter")
    print("=" * 60)
    adata_full = anndata.read_h5ad(INPUT_H5AD)
    sub = adata_full[adata_full.obs["leiden"] == TARGET_CLUSTER].copy()
    n_before = sub.n_obs
    print(f"  CL{TARGET_CLUSTER}: {n_before:,} cells")

    counts = sub.layers["counts"]
    filter_idx = [sub.var_names.get_loc(g) for g in EC_FILTER_GENES if g in sub.var_names]
    print(f"  filter genes: {[g for g in EC_FILTER_GENES if g in sub.var_names]}")
    sub_counts = counts[:, filter_idx].toarray() if sp.issparse(counts) else counts[:, filter_idx]
    is_ec = (sub_counts > 0).any(axis=1)
    sub = sub[is_ec].copy()
    print(
        f"  endothelial (filter > 0): {sub.n_obs:,} / {n_before:,} "
        f"({sub.n_obs / n_before * 100:.1f}%)"
    )
    print(f"  슬라이드: {sub.obs['slide'].value_counts().to_dict()}")
    print(f"  병리: {sub.obs['pathology'].value_counts().to_dict()}")
    print(f"  이전 celltype top:\n{sub.obs['celltype_final'].value_counts().head().to_string()}")

    sub.obs = sub.obs.rename(columns={"leiden": "leiden_parent"})

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
    print("3. EC feature 선택")
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
    marker_out = TBL_DIR / "cl1_marker_genes.csv"
    marker_df.to_csv(marker_out, index=False)
    print(f"  → {marker_out}")

    # ── 프로그램 평균 ──────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("9. 서브클러스터별 프로그램 스코어 평균")
    print("=" * 60)
    score_mat = sub.obs.groupby("leiden_sub", observed=True)[program_keys].mean().loc[sub_cats]
    score_out = TBL_DIR / "cl1_program_scores.csv"
    score_mat.round(4).to_csv(score_out)
    print(score_mat.round(2).to_string())

    # ── 요약 ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("10. 서브클러스터 요약")
    print("=" * 60)
    summary_rows = []
    for cl in sub_cats:
        cells = sub.obs[sub.obs["leiden_sub"] == cl]
        n = len(cells)
        slide_pct = cells["slide"].value_counts(normalize=True) * 100
        path_pct = cells["pathology"].value_counts(normalize=True) * 100
        ct_dist = cells["celltype_final"].value_counts(normalize=True) * 100
        top_markers = marker_df[marker_df["leiden_sub"] == cl].head(8)["gene"].tolist()
        prog_means = score_mat.loc[cl].sort_values(ascending=False)
        top_prog = ", ".join(
            [f"{k.replace('score_', '')}({v:.2f})" for k, v in prog_means.head(3).items()]
        )
        summary_rows.append(
            {
                "leiden_sub": cl,
                "n_cells": n,
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
    summary_out = TBL_DIR / "cl1_summary.csv"
    summary_df.to_csv(summary_out, index=False)
    print(summary_df.to_string(index=False))

    # ── Figures ────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("11. Figures")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    sc.pl.umap(sub, color="leiden_sub", legend_loc="on data", ax=axes[0], show=False, frameon=False)
    axes[0].set_title(f"CL1 sub-Leiden (res={RESOLUTION}, n={n_sub})")
    sc.pl.umap(sub, color="slide", ax=axes[1], show=False, frameon=False)
    axes[1].set_title("Slide")
    sc.pl.umap(sub, color="pathology", ax=axes[2], show=False, frameon=False)
    axes[2].set_title("Pathology")
    plt.tight_layout()
    out = FIG_DIR / "cl1_umap_subclusters.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")

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
    out = FIG_DIR / "cl1_umap_program_scores.png"
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
        title="CL1 서브클러스터 × EC 마커",
        figsize=(max(12, len(feature_genes) * 0.4), max(4, n_sub * 0.45)),
    )
    out = FIG_DIR / "cl1_dotplot_markers.png"
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
    ax.set_title("CL1 서브클러스터 × 프로그램 스코어 (z-score)", fontweight="bold")
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
    out = FIG_DIR / "cl1_program_score_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")

    # ── 저장 ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("12. h5ad 저장")
    print("=" * 60)
    sub.write_h5ad(OUT_H5AD)
    print(f"  → {OUT_H5AD}")
    print("\n완료!")


if __name__ == "__main__":
    main()
