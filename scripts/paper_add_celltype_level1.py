"""final_celltype(60종)을 상위 계통 그룹 celltype_level1(~13종)으로 묶기.

final_celltype 라벨 중 'proliferating'·'emt' 등은 여러 계통(종양/CAF/myeloid)에
걸치므로, subcluster_analysis(계통 출처)를 우선 사용해 정확히 분류한다.
필터 제외(none) 세포는 최초 celltype_final로 계통 부여.

level1 그룹: Tumor / Myoepithelial / T_NK / B_cell / Plasma / Macrophage_Mono /
Dendritic / Mast / Neutrophil / CAF_Fibroblast / Adipocyte / Endothelial / Pericyte

출력:
    04.tables/master_cell_annotation.csv  (celltype_level1 컬럼 추가)
    04.tables/master_celltype_level1_counts.csv
    04.tables/master_level1_mapping.csv
    05.figures/celltype_level1_composition_dcis_mibc.png
    + integrated_qc_passed.h5ad obs에 celltype_level1 기록
"""

from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

# Published verbatim except for the path constants below, which point at this repository
# and at the pool objects the paper's subcluster scripts write. The rules that decide a
# result are the code that produced the paper.

ROOT = Path(__file__).resolve().parent.parent
PROJ = ROOT
MASTER = PROJ / "03.data_processed/subclustered/master_cell_annotation.csv"
INT = PROJ / "03.data_processed/pool_input.h5ad"
TBL = PROJ / "03.data_processed/subclustered"
FIG = PROJ / "03.data_processed/subclustered/figures"

LEVEL1_ORDER = [
    "Tumor",
    "Myoepithelial",
    "T_NK",
    "B_cell",
    "Plasma",
    "Macrophage_Mono",
    "Dendritic",
    "Mast",
    "Neutrophil",
    "CAF_Fibroblast",
    "Adipocyte",
    "Endothelial",
    "Pericyte",
    "Other",
]

DC_TYPES = {"mregDC", "cDC1", "cDC2", "pDC", "langerhans", "moDC", "cDC"}

# 필터 제외(none) 세포: 최초 celltype_final → level1
CTF_TO_L1 = {
    "cancer": "Tumor",
    "myoepithelial cell": "Myoepithelial",
    "CAF": "CAF_Fibroblast",
    "endothelial cell": "Endothelial",
    "B": "B_cell",
    "plasma cell": "Plasma",
    "CD8T": "T_NK",
    "T": "T_NK",
    "Treg": "T_NK",
    "NKT": "T_NK",
    "perivascular NK": "T_NK",
    "M1": "Macrophage_Mono",
    "M2": "Macrophage_Mono",
    "moDC": "Dendritic",
    "cDC": "Dendritic",
    "pDC": "Dendritic",
    "mast cell": "Mast",
    "neutrophil": "Neutrophil",
}


def to_level1(ft, analysis, ctf, qc_status=None):
    f = str(ft).lower()
    a = str(analysis)
    # 오염(contamination) 세포는 분석 계통이 아니라 본래 celltype으로
    if str(qc_status) == "contamination":
        return CTF_TO_L1.get(str(ctf), "Other")
    if a == "CL3+4+10 tumor":
        return "Myoepithelial" if "myoep" in f else "Tumor"
    if a == "CL6+7+13+14 CAF/adipo":
        return "Adipocyte" if ft in ("adipocyte", "preadipo") else "CAF_Fibroblast"
    if a == "CL5+8 B/plasma":
        return "Plasma" if "plasma" in f else "B_cell"
    if a == "CL2+11+12 myeloid/DC":
        return "Dendritic" if ft in DC_TYPES else "Macrophage_Mono"
    if a == "CL0 T/NK":
        return "T_NK"
    if a == "CL9 mast":
        return "Mast"
    if a == "CL1 endothelial":
        return "Pericyte" if "pericyte" in f else "Endothelial"
    # none (필터 제외) → 최초 celltype_final
    return CTF_TO_L1.get(str(ctf), "Other")


def main():
    m = pd.read_csv(MASTER)
    m["celltype_level1"] = [
        to_level1(ft, an, ctf, qs)
        for ft, an, ctf, qs in zip(
            m["final_celltype"],
            m["subcluster_analysis"],
            m["celltype_final"],
            m["qc_status"],
            strict=True,
        )
    ]
    m.to_csv(MASTER, index=False)
    print(f"master_cell_annotation.csv에 celltype_level1 추가 ({len(m):,} cells)")

    # ── level1 counts ──────────────────────────────────────────────────────
    cnt = m["celltype_level1"].value_counts().reindex(LEVEL1_ORDER).dropna().astype(int)
    cdf = cnt.reset_index()
    cdf.columns = ["celltype_level1", "n_cells"]
    cdf["pct"] = (cdf["n_cells"] / len(m) * 100).round(2)
    cdf.to_csv(TBL / "master_celltype_level1_counts.csv", index=False)
    print("\n=== celltype_level1 분포 ===")
    print(cdf.to_string(index=False))

    # ── level1 → final_celltype 매핑 ───────────────────────────────────────
    mp = (
        m.groupby("celltype_level1")["final_celltype"]
        .agg(lambda x: ", ".join(sorted(x.unique())))
        .reset_index()
    )
    mp["n_final_celltype"] = m.groupby("celltype_level1")["final_celltype"].nunique().values
    mp = mp[["celltype_level1", "n_final_celltype", "final_celltype"]]
    mp.to_csv(TBL / "master_level1_mapping.csv", index=False)
    print("\n=== level1 → final_celltype 매핑 ===")
    for _, r in mp.iterrows():
        print(
            f"  {r['celltype_level1']:<16} ({r['n_final_celltype']:>2}): {r['final_celltype'][:110]}"
        )

    # ── 그림: level1 구성 DCIS vs MIBC (코어 평균 비율) ─────────────────────
    m["core"] = m["slide"] + "_" + m["core_id"].astype(str)
    ct = (
        m.groupby(["core", "pathology", "celltype_level1"], observed=True)
        .size()
        .reset_index(name="n")
    )
    ct["total"] = ct.groupby("core")["n"].transform("sum")
    ct["prop"] = ct["n"] / ct["total"] * 100
    order = [g for g in LEVEL1_ORDER if g in m["celltype_level1"].unique()]
    means = (
        ct.groupby(["pathology", "celltype_level1"], observed=True)["prop"]
        .mean()
        .unstack(fill_value=0)
        .reindex(columns=order)
    )
    colors = plt.cm.tab20(np.linspace(0, 1, len(order)))
    fig, ax = plt.subplots(figsize=(7, 7))
    bottom = np.zeros(2)
    xlab = ["DCIS", "MIBC"]
    for i, g in enumerate(order):
        vals = means.loc[["dcis", "mibc"], g].values
        ax.bar(xlab, vals, bottom=bottom, color=colors[i], label=g)
        for xi, (v, b) in enumerate(zip(vals, bottom, strict=True)):
            if v >= 3:
                ax.text(xi, b + v / 2, f"{v:.0f}", ha="center", va="center", fontsize=7)
        bottom += vals
    ax.set_ylabel("코어 평균 구성 비율 (%)")
    ax.set_title("celltype_level1 구성: DCIS vs MIBC (코어 평균)")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "celltype_level1_composition_dcis_mibc.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("\n저장: celltype_level1_composition_dcis_mibc.png")

    # ── integrated h5ad obs에 기록 ─────────────────────────────────────────
    a = ad.read_h5ad(INT)
    l1 = [
        to_level1(ft, an, ctf, qs)
        for ft, an, ctf, qs in zip(
            a.obs["final_celltype"].astype(str),
            a.obs["subcluster_analysis"].astype(str),
            a.obs["celltype_final"].astype(str),
            a.obs["qc_status"].astype(str),
            strict=True,
        )
    ]
    a.obs["celltype_level1"] = pd.Categorical(l1)
    a.write_h5ad(INT)
    print("integrated_qc_passed.h5ad obs에 celltype_level1 기록 완료.\n완료!")


if __name__ == "__main__":
    main()
