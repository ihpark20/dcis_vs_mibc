"""전체 서브클러스터를 annotation 기반으로 병합 (재사용 모듈).

각 서브클러스터링 결과(summary CSV + h5ad)에 대해 3가지 후처리:
  1. fragment 플래그: n_cells < max(200, 0.5% of pool)
  2. contamination 플래그: dominant_prev_celltype가 해당 풀의 정당 계통(target)이
     아니고 비율 >= CONTAM_PCT
  3. 병합: leiden_sub → subtype 로 묶음 (subtype_merged), qc_status 동반

각 h5ad obs에 subtype_merged / qc_status 기록하고, 전체 매핑·요약 테이블 + 그림 생성.

출력:
    04.tables/subcluster_merged_mapping.csv   (모든 leiden_sub × 병합 라벨)
    04.tables/subcluster_merged_summary.csv   (dataset × subtype_merged)
    05.figures/subcluster_merge_overview.png
    + 각 h5ad obs에 subtype_merged, qc_status 추가
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

# Published verbatim except for the path constants below and the dataset list, which point
# at this repository and at the seven pools the paper's annotation used. The rules that
# decide a result — the fragment threshold, the contamination rule and the merge — are the
# code that produced the paper.

ROOT = Path(__file__).resolve().parent.parent
PROJ = ROOT
TBL = PROJ / "03.data_processed/subclustered"
H5 = PROJ / "03.data_processed/subclustered"
FIG = PROJ / "03.data_processed/subclustered/figures"

FRAG_ABS = 200
FRAG_REL = 0.005
CONTAM_PCT = 40.0

# dataset config: summary CSV, h5ad, 정당 계통(celltype_final), 표시 이름
DATASETS = [
    {
        "name": "CL1 endothelial",
        "summary": "cl1_subtype_annotation.csv",
        "h5ad": "cl1_endothelial_subclustered.h5ad",
        "target": {"endothelial cell"},
    },
    {
        "name": "CL3+4+10 tumor",
        "summary": "cl3_4_10_combined_summary.csv",
        "h5ad": "cl3_4_10_combined_subclustered.h5ad",
        "target": {"cancer", "myoepithelial cell"},
    },
    {
        "name": "CL6+7+13+14 CAF+adipo",
        "summary": "cl6_7_13_14_caf_adipo_summary.csv",
        "h5ad": "cl6_7_13_14_caf_adipo_subclustered.h5ad",
        "target": {"CAF"},
    },
    {
        "name": "CL5+8 B/plasma",
        "summary": "cl5_8_bplasma_summary.csv",
        "h5ad": "cl5_8_bplasma_subclustered.h5ad",
        "target": {"B", "plasma cell"},
    },
    {
        "name": "CL9 mast",
        "summary": "cl9_mast_summary.csv",
        "h5ad": "cl9_mast_subclustered.h5ad",
        "target": {"mast cell"},
    },
    {
        "name": "CL2+11+12 myeloid/DC",
        "summary": "cl2_11_12_myeloid_dc_summary.csv",
        "h5ad": "cl2_11_12_myeloid_dc_subclustered.h5ad",
        "target": {"M1", "M2", "moDC", "cDC", "pDC"},
    },
]


def classify(summary, target, pool_total):
    """leiden_sub별 qc_status + subtype_merged 결정."""
    frag_thr = max(FRAG_ABS, FRAG_REL * pool_total)
    out = []
    for _, r in summary.iterrows():
        sub = str(r["leiden_sub"])
        n = int(r["n_cells"])
        subtype = r["subtype"]
        dom = r.get("dominant_prev_celltype", None)
        dom_pct = r.get("dominant_prev_pct", 0.0)
        is_frag = (
            bool(r["is_fragment"])
            if "is_fragment" in r and not pd.isna(r["is_fragment"])
            else (n < frag_thr)
        )
        is_contam = isinstance(dom, str) and dom not in target and float(dom_pct or 0) >= CONTAM_PCT
        if is_contam:
            status = "contamination"
            merged = f"contamination:{dom}"
        elif is_frag:
            status = "fragment"
            merged = subtype  # fragment도 생물학적으로는 subtype에 귀속
        else:
            status = "clean"
            merged = subtype
        out.append(
            {
                "leiden_sub": sub,
                "n_cells": n,
                "subtype": subtype,
                "dominant_prev_celltype": dom,
                "dominant_prev_pct": dom_pct,
                "qc_status": status,
                "subtype_merged": merged,
            }
        )
    return pd.DataFrame(out)


def main():
    all_map = []
    summary_rows = []
    for ds in DATASETS:
        s = pd.read_csv(TBL / ds["summary"])
        s["leiden_sub"] = s["leiden_sub"].astype(str)
        pool_total = int(s["n_cells"].sum())
        cls = classify(s, ds["target"], pool_total)
        cls.insert(0, "dataset", ds["name"])
        all_map.append(cls)

        n_sub = len(cls)
        n_merged = cls.loc[cls["qc_status"] != "contamination", "subtype_merged"].nunique()
        n_clean = (cls["qc_status"] == "clean").sum()
        n_frag = (cls["qc_status"] == "fragment").sum()
        n_contam = (cls["qc_status"] == "contamination").sum()
        cells_contam = int(cls.loc[cls["qc_status"] == "contamination", "n_cells"].sum())
        print(
            f"{ds['name']:<24} {n_sub:>3} sub → {n_merged:>2} subtype "
            f"| clean {n_clean}, fragment {n_frag}, contam {n_contam} "
            f"(오염 {cells_contam:,}/{pool_total:,} = {cells_contam / pool_total * 100:.1f}%)"
        )
        summary_rows.append(
            {
                "dataset": ds["name"],
                "n_subclusters": n_sub,
                "pool_cells": pool_total,
                "n_merged_subtypes": n_merged,
                "n_clean": n_clean,
                "n_fragment": n_frag,
                "n_contamination": n_contam,
                "pct_contamination_cells": round(cells_contam / pool_total * 100, 1),
            }
        )

        # ── h5ad obs에 기록 ──
        h5path = H5 / ds["h5ad"]
        if h5path.exists():
            a = ad.read_h5ad(h5path)
            key = "leiden_sub" if "leiden_sub" in a.obs.columns else None
            if key:
                m_merged = dict(zip(cls["leiden_sub"], cls["subtype_merged"], strict=True))
                m_status = dict(zip(cls["leiden_sub"], cls["qc_status"], strict=True))
                a.obs["subtype_merged"] = a.obs[key].astype(str).map(m_merged).astype("category")
                a.obs["qc_status"] = a.obs[key].astype(str).map(m_status).astype("category")
                a.write_h5ad(h5path)

    mapping = pd.concat(all_map, ignore_index=True)
    mapping.to_csv(TBL / "subcluster_merged_mapping.csv", index=False)
    msum = pd.DataFrame(summary_rows)
    msum.to_csv(TBL / "subcluster_merged_summary.csv", index=False)
    print("\n저장: subcluster_merged_mapping.csv, subcluster_merged_summary.csv")
    print("\n=== dataset별 병합 요약 ===")
    print(msum.to_string(index=False))

    # ── 그림: before(서브클러스터 수) vs after(subtype 수) + 오염% ──────────
    msum = msum.sort_values("n_subclusters", ascending=True)
    y = np.arange(len(msum))
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    ax = axes[0]
    ax.barh(y - 0.2, msum["n_subclusters"], 0.4, color="#bdbdbd", label="병합 전 (leiden_sub)")
    ax.barh(y + 0.2, msum["n_merged_subtypes"], 0.4, color="#3182bd", label="병합 후 (subtype)")
    for yi, r in zip(y, msum.itertuples(), strict=True):
        ax.text(r.n_subclusters + 0.5, yi - 0.2, str(r.n_subclusters), va="center", fontsize=8)
        ax.text(
            r.n_merged_subtypes + 0.5,
            yi + 0.2,
            str(r.n_merged_subtypes),
            va="center",
            fontsize=8,
            color="#08519c",
        )
    ax.set_yticks(y)
    ax.set_yticklabels(msum["dataset"], fontsize=9)
    ax.set_xlabel("클러스터 수")
    ax.set_title("(A) annotation 병합: leiden_sub → subtype", fontsize=12, fontweight="bold")
    ax.legend(loc="lower right")

    ax = axes[1]
    cl = msum.set_index("dataset")[["n_clean", "n_fragment", "n_contamination"]].loc[
        msum["dataset"]
    ]
    left = np.zeros(len(msum))
    for col, color, lab in [
        ("n_clean", "#74c476", "clean"),
        ("n_fragment", "#fdae6b", "fragment"),
        ("n_contamination", "#d62728", "contamination"),
    ]:
        ax.barh(y, cl[col].values, left=left, color=color, label=lab, height=0.6)
        left += cl[col].values
    ax.set_yticks(y)
    ax.set_yticklabels(msum["dataset"], fontsize=9)
    ax.set_xlabel("서브클러스터 수")
    ax.set_title("(B) 서브클러스터 QC 상태", fontsize=12, fontweight="bold")
    ax.legend(loc="lower right")
    fig.suptitle("서브클러스터 annotation 기반 병합 개요", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "subcluster_merge_overview.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("저장: subcluster_merge_overview.png")
    print("\n완료! (각 h5ad obs에 subtype_merged, qc_status 추가됨)")


if __name__ == "__main__":
    main()
