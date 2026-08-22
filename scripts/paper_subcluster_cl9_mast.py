"""CL9 Mast cell 서브클러스터링.

CL9 (mast cell 우세, 8,559 — CAF/neutrophil 오염 다수)를 mast 마커 하나라도
양성인 세포만 필터 후 재클러스터링. subtype/마커/프로그램 + 코어 dominance
+ 서브클러스터별 DCIS vs MIBC 비교.

전략:
  1. CL9 subset (QC FAIL 제외)
  2. mast 마커(>0) 하나라도 양성 필터
  3. mast focused feature로 PCA→Harmony(slide)→Leiden
  4. 프로그램 스코어 + 코어 dominance(Top1/2/3) + DCIS vs MIBC(코어 단위 검정)

출력:
    03.data_processed/cl9_mast_subclustered.h5ad
    04.tables/cl9_mast_marker_genes.csv
    04.tables/cl9_mast_program_scores.csv
    04.tables/cl9_mast_summary.csv             (subtype/마커/프로그램/dominance/DCIS·MIBC%)
    04.tables/cl9_mast_dcis_vs_mibc_test.csv
    05.figures/cl9_mast_umap_subclusters.png
    05.figures/cl9_mast_umap_program_scores.png
    05.figures/cl9_mast_dotplot_markers.png
    05.figures/cl9_mast_program_score_heatmap.png
    05.figures/cl9_mast_dcis_vs_mibc_boxplot.png
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
from scipy import stats  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

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
OUT_H5AD = PROJ / "03.data_processed/subclustered/cl9_mast_subclustered.h5ad"
FIG_DIR = PROJ / "03.data_processed/subclustered/figures"
TBL_DIR = PROJ / "03.data_processed/subclustered"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(FIG_DIR)
sc.settings.verbosity = 1

TARGET_CLUSTERS = ["9"]
RESOLUTION = 0.3
RANDOM_STATE = 42
N_PCS = 12
LABEL = "cl9_mast"
QC_FAIL = {("TMA1", "D4"), ("TMA1", "D8")}
PATH_COLORS = {"dcis": "#4ECDC4", "mibc": "#FF6B6B"}

# mast 마커 하나라도 양성 → 유지
FILTER_GENES = ["TPSAB1", "CPA3", "KIT", "FCER1A", "GATA2", "CTSG"]

FEATURE_GENES = [
    # mast core
    "TPSAB1",
    "CPA3",
    "KIT",
    "FCER1A",
    "GATA2",
    "CTSG",
    "CD9",
    # effector / activation
    "IL4",
    "IL13",
    "TNF",
    "CD69",
    "IL6",
    "CCL5",
    "PTGS2",
    "AREG",
    "CSF1",
    # IFN response
    "STAT1",
    "MX1",
    "IFIT3",
    "IRF1",
    "ISG15",
    # proliferation
    "MKI67",
    "TOP2A",
    "CENPF",
    "HMGA1",
    # contamination 식별용 (CAF/neutrophil/myeloid)
    "PDGFRB",
    "LUM",
    "FN1",
    "LYZ",
    "CD68",
    "S100A8",
    "FCGR3B",
]

PROGRAM_MARKERS = {
    "mast_core": ["TPSAB1", "CPA3", "KIT", "FCER1A", "GATA2", "CTSG", "CD9"],
    "effector_activation": ["IL4", "IL13", "TNF", "CD69", "IL6", "CCL5", "PTGS2", "AREG", "CSF1"],
    "ifn_response": ["STAT1", "MX1", "IFIT3", "IRF1", "ISG15"],
    "proliferating": ["MKI67", "TOP2A", "CENPF", "HMGA1"],
}
SUBTYPE_PROGS = ["mast_core", "effector_activation", "ifn_response", "proliferating"]


def main():
    print("1. CL9 subset + filter")
    adata_full = anndata.read_h5ad(INPUT_H5AD)
    sub = adata_full[adata_full.obs["leiden"].isin(TARGET_CLUSTERS)].copy()
    for s, c in QC_FAIL:
        sub = sub[~((sub.obs["slide"] == s) & (sub.obs["core_id"] == c))].copy()
    n_before = sub.n_obs

    counts = sub.layers["counts"]
    filter_idx = [sub.var_names.get_loc(g) for g in FILTER_GENES if g in sub.var_names]
    print(f"  filter genes: {[g for g in FILTER_GENES if g in sub.var_names]}")
    sub_counts = counts[:, filter_idx].toarray() if sp.issparse(counts) else counts[:, filter_idx]
    is_pass = (sub_counts > 0).any(axis=1)
    sub = sub[is_pass].copy()
    print(f"  filter > 0: {sub.n_obs:,} / {n_before:,} ({sub.n_obs / n_before * 100:.1f}%)")
    print(f"  병리: {sub.obs['pathology'].value_counts().to_dict()}")
    print(f"  이전 celltype: {sub.obs['celltype_final'].value_counts().head(5).to_dict()}")

    sub.obs = sub.obs.drop(columns=["leiden"])

    print("\n2. 정규화")
    sub.X = sub.layers["counts"].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sub.layers["lognorm"] = sub.X.copy()

    feature_genes = [g for g in FEATURE_GENES if g in sub.var_names]
    print(f"3. Feature: {len(feature_genes)}/{len(FEATURE_GENES)}")
    sub.var["highly_variable"] = sub.var_names.isin(feature_genes)

    print("4. Scale + PCA")
    sc.pp.scale(sub, max_value=10)
    n_comps = min(N_PCS, len(feature_genes) - 1)
    sc.tl.pca(sub, n_comps=n_comps, use_highly_variable=True)

    print("5. Harmony")
    harmony_out = hm.run_harmony(
        sub.obsm["X_pca"].copy(),
        sub.obs[["slide"]].copy(),
        vars_use=["slide"],
        random_state=RANDOM_STATE,
        max_iter_harmony=30,
    )
    sub.obsm["X_pca_harmony"] = np.array(harmony_out.Z_corr)

    print(f"6. Neighbors + UMAP + Leiden (res={RESOLUTION})")
    sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=n_comps)
    sc.tl.umap(sub, random_state=RANDOM_STATE)
    sc.tl.leiden(sub, resolution=RESOLUTION, random_state=RANDOM_STATE, key_added="leiden_sub")
    n_sub = sub.obs["leiden_sub"].nunique()
    print(f"  → {n_sub} 서브클러스터")
    print(sub.obs["leiden_sub"].value_counts().sort_index().to_string())

    print("7. 프로그램 스코어")
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

    print("8. 마커 (Wilcoxon)")
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
        u = sub_score.uns["rank_sub"]
        for r, (g, s, p, lfc) in enumerate(
            zip(
                u["names"][cl],
                u["scores"][cl],
                u["pvals_adj"][cl],
                u["logfoldchanges"][cl],
                strict=False,
            ),
            start=1,
        ):
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

    score_mat = sub.obs.groupby("leiden_sub", observed=True)[program_keys].mean().loc[sub_cats]
    score_mat.round(4).to_csv(TBL_DIR / f"{LABEL}_program_scores.csv")

    print("9. 요약 + dominance + DCIS/MIBC%")
    sub.obs["core"] = sub.obs["slide"].astype(str) + "_" + sub.obs["core_id"].astype(str)
    core_path = sub.obs.drop_duplicates("core").set_index("core")["pathology"].to_dict()
    rows = []
    for cl in sub_cats:
        cells = sub.obs[sub.obs["leiden_sub"] == cl]
        n = len(cells)
        path_pct = cells["pathology"].value_counts(normalize=True) * 100
        ct_dist = cells["celltype_final"].value_counts(normalize=True) * 100
        top_markers = marker_df[marker_df["leiden_sub"] == cl].head(8)["gene"].tolist()
        prog_means = score_mat.loc[cl].sort_values(ascending=False)
        top_prog = ", ".join(
            [f"{k.replace('score_', '')}({v:.2f})" for k, v in prog_means.head(3).items()]
        )
        subtype_label = max(
            {p: score_mat.loc[cl].get(f"score_{p}", -np.inf) for p in SUBTYPE_PROGS},
            key=lambda p: score_mat.loc[cl].get(f"score_{p}", -np.inf),
        )
        cc = cells["core"].value_counts()
        share = (cc / cc.sum() * 100).values
        cum = np.cumsum(share)
        tops = cc.index.tolist()
        hhi = float(np.sum((share / 100) ** 2))
        rows.append(
            {
                "leiden_sub": cl,
                "n_cells": n,
                "subtype": subtype_label,
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
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(TBL_DIR / f"{LABEL}_summary.csv", index=False)
    print(
        summary_df[
            [
                "leiden_sub",
                "n_cells",
                "subtype",
                "pct_mibc",
                "pct_dcis",
                "dominant_prev_celltype",
                "top1_pct",
                "top3_cum_pct",
                "eff_cores",
                "top8_markers",
            ]
        ].to_string(index=False)
    )

    # ── 서브클러스터별 DCIS vs MIBC (코어 단위 비율 검정) ───────────────────
    print("\n10. DCIS vs MIBC 비교 (코어 단위 비율, Mann-Whitney + BH FDR)")
    ct = (
        sub.obs.groupby(["core", "pathology", "leiden_sub"], observed=True)
        .size()
        .reset_index(name="n")
    )
    ct["total"] = ct.groupby("core")["n"].transform("sum")
    ct["prop"] = ct["n"] / ct["total"]
    prop = ct.pivot_table(
        index=["core", "pathology"],
        columns="leiden_sub",
        values="prop",
        fill_value=0,
        observed=True,
    ).reset_index()
    prop.columns.name = None
    dcis = prop[prop["pathology"] == "dcis"]
    mibc = prop[prop["pathology"] == "mibc"]
    test_rows = []
    for cl in sub_cats:
        if cl not in prop.columns:
            continue
        u, p = stats.mannwhitneyu(dcis[cl], mibc[cl], alternative="two-sided")
        test_rows.append(
            {
                "leiden_sub": cl,
                "subtype": summary_df.set_index("leiden_sub").loc[cl, "subtype"],
                "dcis_mean_prop": round(dcis[cl].mean(), 4),
                "mibc_mean_prop": round(mibc[cl].mean(), 4),
                "log2fc_mibc_dcis": round(
                    np.log2((mibc[cl].mean() + 1e-4) / (dcis[cl].mean() + 1e-4)), 2
                ),
                "p_value": p,
            }
        )
    test = pd.DataFrame(test_rows)
    test["p_adj"] = multipletests(test["p_value"], method="fdr_bh")[1]
    test["sig"] = ["*" if x < 0.05 else "ns" for x in test["p_adj"]]
    test = test.sort_values("p_adj")
    test.to_csv(TBL_DIR / f"{LABEL}_dcis_vs_mibc_test.csv", index=False)
    print(test.to_string(index=False))

    print("\n11. Figures")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    sc.pl.umap(sub, color="leiden_sub", legend_loc="on data", ax=axes[0], show=False, frameon=False)
    axes[0].set_title(f"CL9 mast sub (res={RESOLUTION}, n={n_sub})")
    sc.pl.umap(sub, color="celltype_final", ax=axes[1], show=False, frameon=False)
    axes[1].set_title("이전 celltype")
    sc.pl.umap(sub, color="pathology", ax=axes[2], show=False, frameon=False)
    axes[2].set_title("Pathology")
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{LABEL}_umap_subclusters.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    n_col = min(4, len(program_keys))
    n_row = int(np.ceil(len(program_keys) / n_col))
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
        figsize=(max(12, len(feature_genes) * 0.4), max(4, n_sub * 0.5)),
    )
    fig_dot.savefig(FIG_DIR / f"{LABEL}_dotplot_markers.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    score_z = ((score_mat - score_mat.mean()) / score_mat.std()).clip(-3, 3)
    fig, ax = plt.subplots(figsize=(max(7, len(program_keys) * 0.8), max(4, n_sub * 0.45)))
    im = ax.imshow(score_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(program_keys)))
    ax.set_xticklabels([k.replace("score_", "") for k in program_keys], rotation=45, ha="right")
    ax.set_yticks(range(n_sub))
    ax.set_yticklabels([f"sub{c}" for c in sub_cats])
    plt.colorbar(im, ax=ax, label="z", shrink=0.7)
    for i in range(n_sub):
        for j in range(len(program_keys)):
            ax.text(j, i, f"{score_mat.values[i, j]:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("CL9 mast × 프로그램 (z)", fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG_DIR / f"{LABEL}_program_score_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # DCIS vs MIBC boxplot
    x = np.arange(len(sub_cats))
    w = 0.38
    fig, ax = plt.subplots(figsize=(max(10, len(sub_cats) * 1.1), 6))
    bp_d = ax.boxplot(
        [dcis[cl].values for cl in sub_cats if cl in prop.columns],
        positions=x - w / 2,
        widths=w * 0.85,
        patch_artist=True,
        boxprops=dict(facecolor=PATH_COLORS["dcis"], alpha=0.7),
        medianprops=dict(color="black"),
        showfliers=False,
    )
    bp_m = ax.boxplot(
        [mibc[cl].values for cl in sub_cats if cl in prop.columns],
        positions=x + w / 2,
        widths=w * 0.85,
        patch_artist=True,
        boxprops=dict(facecolor=PATH_COLORS["mibc"], alpha=0.7),
        medianprops=dict(color="black"),
        showfliers=False,
    )
    rng = np.random.default_rng(0)
    tmap = test.set_index("leiden_sub")
    for i, cl in enumerate(sub_cats):
        if cl not in prop.columns:
            continue
        ax.scatter(
            rng.normal(x[i] - w / 2, 0.05, len(dcis)),
            dcis[cl],
            s=10,
            color="#2BA5A0",
            alpha=0.5,
            zorder=3,
        )
        ax.scatter(
            rng.normal(x[i] + w / 2, 0.05, len(mibc)),
            mibc[cl],
            s=10,
            color="#E04545",
            alpha=0.5,
            zorder=3,
        )
        ymax = max(dcis[cl].max(), mibc[cl].max())
        ax.text(
            x[i],
            ymax + 0.02,
            tmap.loc[cl, "sig"],
            ha="center",
            fontsize=9,
            color="red" if tmap.loc[cl, "sig"] != "ns" else "gray",
            fontweight="bold",
        )
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"sub{cl}\n{summary_df.set_index('leiden_sub').loc[cl, 'subtype']}" for cl in sub_cats],
        fontsize=8,
    )
    ax.set_ylabel("코어 내 비율 (mast 중)")
    ax.set_title("CL9 mast 서브클러스터: DCIS vs MIBC (코어 단위, BH-FDR)")
    ax.legend([bp_d["boxes"][0], bp_m["boxes"][0]], ["DCIS", "MIBC"], loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"{LABEL}_dcis_vs_mibc_boxplot.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print("\n12. h5ad 저장")
    sub.write_h5ad(OUT_H5AD)
    print(f"  → {OUT_H5AD}\n완료!")


if __name__ == "__main__":
    main()
