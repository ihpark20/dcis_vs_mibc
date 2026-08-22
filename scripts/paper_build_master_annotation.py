"""전체 데이터셋 master cell annotation 테이블 생성.

각 최초 Leiden 클러스터(CL0–CL14)를 하나의 (통합 우선) 서브클러스터링 결과에
매핑하고, 세포 단위로 다음 정보를 통합한다:
  - orig_cluster        : 최초 Leiden 클러스터 (CL0–CL14)
  - subcluster_analysis : 사용한 서브클러스터링 분석 (통합 우선)
  - covered_clusters    : 그 분석이 포함한 최초 클러스터 세트
  - leiden_sub          : 그 분석에서의 서브클러스터 라벨
  - subtype_merged      : annotation 기반 병합 라벨 (merge_subclusters)
  - qc_status           : clean / fragment / contamination
  - celltype_final      : 최초 celltype 주석 (broad)
  - final_celltype      : 최종 주석 (clean/fragment면 subtype_merged, 아니면 celltype_final)

최초 클러스터가 여러 서브클러스터링에 포함된 경우 통합 버전을 채택:
  CL2 → CL2+11+12, CL3/4/10 → 통합, CL6/7/13/14 → 통합, CL5/8 → 통합, CL11/12 → CL2+11+12.

출력:
    04.tables/master_cell_annotation.csv
    04.tables/master_annotation_legend.csv
    04.tables/master_final_celltype_counts.csv
    + integrated_qc_passed.h5ad obs에 master 컬럼 기록
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core_exclusion import EXCLUDED_CORES, keep_mask  # noqa: E402

# Published verbatim except for the path constants below, which point at this repository
# and at the pool objects the paper's subcluster scripts write. The rules that decide a
# result are the code that produced the paper.

ROOT = Path(__file__).resolve().parent.parent
PROJ = ROOT
H5 = PROJ / "03.data_processed/subclustered"
TBL = PROJ / "03.data_processed/subclustered"
INT = PROJ / "03.data_processed/pool_input.h5ad"

FRAG_ABS = 200
FRAG_REL = 0.005
CONTAM_PCT = 40.0

# 채택한 서브클러스터링 (통합 우선). covers = 포함 최초 클러스터.
CHOSEN = [
    {
        "name": "CL0 T/NK",
        "h5ad": "cl0_tcell_only_subclustered.h5ad",
        "covers": ["CL0"],
        "summary": "cl0_tcell_only_summary.csv",
        "target": {"CD8T", "T", "Treg", "NKT", "perivascular NK"},
    },
    {"name": "CL1 endothelial", "h5ad": "cl1_endothelial_subclustered.h5ad", "covers": ["CL1"]},
    {
        "name": "CL2+11+12 myeloid/DC",
        "h5ad": "cl2_11_12_myeloid_dc_subclustered.h5ad",
        "covers": ["CL2", "CL11", "CL12"],
    },
    {
        "name": "CL3+4+10 tumor",
        "h5ad": "cl3_4_10_combined_subclustered.h5ad",
        "covers": ["CL3", "CL4", "CL10"],
    },
    {"name": "CL5+8 B/plasma", "h5ad": "cl5_8_bplasma_subclustered.h5ad", "covers": ["CL5", "CL8"]},
    {
        "name": "CL6+7+13+14 CAF/adipo",
        "h5ad": "cl6_7_13_14_caf_adipo_subclustered.h5ad",
        "covers": ["CL6", "CL7", "CL13", "CL14"],
    },
    {"name": "CL9 mast", "h5ad": "cl9_mast_subclustered.h5ad", "covers": ["CL9"]},
]


def derive_cl0(summary_path, target):
    """subtype_merged 없는 CL0: top_programs에서 subtype 유도 + fragment/contam 플래그."""
    s = pd.read_csv(TBL / summary_path)
    s["leiden_sub"] = s["leiden_sub"].astype(str)
    pool = int(s["n_cells"].sum())
    frag_thr = max(FRAG_ABS, FRAG_REL * pool)
    sub_map, status_map = {}, {}
    for _, r in s.iterrows():
        subtype = str(r["top_programs"]).split("(")[0].strip()
        dom = r.get("dominant_prev_celltype")
        dom_pct = float(r.get("dominant_prev_pct", 0) or 0)
        if isinstance(dom, str) and dom not in target and dom_pct >= CONTAM_PCT:
            status_map[r["leiden_sub"]] = "contamination"
            sub_map[r["leiden_sub"]] = f"contamination:{dom}"
        else:
            status_map[r["leiden_sub"]] = "fragment" if r["n_cells"] < frag_thr else "clean"
            sub_map[r["leiden_sub"]] = subtype
    return sub_map, status_map


def main():
    print("integrated 로드...")
    a = ad.read_h5ad(INT)
    obs = a.obs
    master = pd.DataFrame(index=obs.index)
    master["orig_cluster"] = "CL" + obs["leiden"].astype(str)
    master["celltype_final"] = obs["celltype_final"].astype(str)
    master["celltype_group"] = obs["celltype_group"].astype(str)
    master["subcluster_analysis"] = "none"
    master["covered_clusters"] = ""
    master["leiden_sub"] = pd.NA
    master["subtype_merged"] = pd.NA
    master["qc_status"] = pd.NA

    for ds in CHOSEN:
        h5 = H5 / ds["h5ad"]
        sub = ad.read_h5ad(h5, backed="r")
        sobs = sub.obs
        idx = master.index.intersection(sobs.index)
        ls = sobs.loc[idx, "leiden_sub"].astype(str)
        if "subtype_merged" in sobs.columns:
            sm = sobs.loc[idx, "subtype_merged"].astype(str)
            qs = sobs.loc[idx, "qc_status"].astype(str)
        else:  # CL0
            sub_map, status_map = derive_cl0(ds["summary"], ds["target"])
            sm = ls.map(sub_map)
            qs = ls.map(status_map)
        master.loc[idx, "subcluster_analysis"] = ds["name"]
        master.loc[idx, "covered_clusters"] = "+".join(ds["covers"])
        master.loc[idx, "leiden_sub"] = ls.values
        master.loc[idx, "subtype_merged"] = sm.values
        master.loc[idx, "qc_status"] = qs.values
        print(f"  {ds['name']:<24} 매핑 {len(idx):>7,} cells")

    # final_celltype: clean/fragment면 subtype_merged, 아니면 celltype_final(broad)
    use_sub = master["qc_status"].isin(["clean", "fragment"]) & master["subtype_merged"].notna()
    master["final_celltype"] = np.where(use_sub, master["subtype_merged"], master["celltype_final"])

    # 좌표·메타
    master.insert(0, "cell_id", master.index)
    for c in ["slide", "core_id", "pathology", "x_centroid", "y_centroid"]:
        master[c] = obs[c].values
    cols = [
        "cell_id",
        "slide",
        "core_id",
        "pathology",
        "x_centroid",
        "y_centroid",
        "orig_cluster",
        "celltype_final",
        "celltype_group",
        "subcluster_analysis",
        "covered_clusters",
        "leiden_sub",
        "subtype_merged",
        "qc_status",
        "final_celltype",
    ]
    master_full = master[cols + []].copy()
    master_full["analysis_exclude"] = ~keep_mask(master_full).values
    # 분석 제외 코어(D4/D8/D6) 제거한 clean master 저장
    master = master_full[~master_full["analysis_exclude"]][cols]
    master.to_csv(TBL / "master_cell_annotation.csv", index=False)
    n_excl = int(master_full["analysis_exclude"].sum())
    print(
        f"\n저장: master_cell_annotation.csv ({len(master):,} cells; "
        f"분석 제외 코어 {sorted({(s, c) for (s, c) in EXCLUDED_CORES})} → {n_excl:,} cells 제외)"
    )

    # ── legend: 분석 → 포함 클러스터, 세포수, subtype 수 ────────────────────
    leg = []
    for ds in CHOSEN:
        m = master[master["subcluster_analysis"] == ds["name"]]
        n_subtype = m.loc[m["qc_status"] != "contamination", "subtype_merged"].nunique()
        leg.append(
            {
                "subcluster_analysis": ds["name"],
                "covered_clusters": "+".join(ds["covers"]),
                "n_cells": len(m),
                "n_leiden_sub": m["leiden_sub"].nunique(),
                "n_subtype_merged": n_subtype,
                "n_clean": int((m["qc_status"] == "clean").sum()) if len(m) else 0,
                "pct_contamination": round((m["qc_status"] == "contamination").mean() * 100, 1)
                if len(m)
                else 0,
            }
        )
    legend = pd.DataFrame(leg)
    n_none = int((master["subcluster_analysis"] == "none").sum())
    legend.loc[len(legend)] = {
        "subcluster_analysis": "(none: 필터 제외)",
        "covered_clusters": "-",
        "n_cells": n_none,
        "n_leiden_sub": 0,
        "n_subtype_merged": 0,
        "n_clean": 0,
        "pct_contamination": 0,
    }
    legend.to_csv(TBL / "master_annotation_legend.csv", index=False)
    print("\n=== 분석 legend ===")
    print(legend.to_string(index=False))

    # ── final_celltype 분포 ────────────────────────────────────────────────
    fc = master["final_celltype"].value_counts().reset_index()
    fc.columns = ["final_celltype", "n_cells"]
    fc["pct"] = (fc["n_cells"] / len(master) * 100).round(2)
    fc.to_csv(TBL / "master_final_celltype_counts.csv", index=False)
    print(f"\n=== final_celltype 종류: {len(fc)} ===")
    print(fc.head(40).to_string(index=False))

    # ── integrated h5ad obs에 기록 (전체 세포; 제외 코어는 analysis_exclude 플래그) ──
    mf = master_full.reindex(a.obs.index)
    for c in [
        "orig_cluster",
        "subcluster_analysis",
        "covered_clusters",
        "leiden_sub",
        "subtype_merged",
        "qc_status",
        "final_celltype",
    ]:
        vals = mf[c].astype("object").where(mf[c].notna(), "NA").astype(str)
        a.obs[f"master_{c}" if c == "leiden_sub" else c] = pd.Categorical(vals.values)
    a.obs["analysis_exclude"] = pd.Categorical(mf["analysis_exclude"].astype(str).values)
    a.write_h5ad(INT)
    print("\nintegrated_qc_passed.h5ad obs에 master 컬럼 + analysis_exclude 기록 완료.")
    print("완료!")


if __name__ == "__main__":
    main()
