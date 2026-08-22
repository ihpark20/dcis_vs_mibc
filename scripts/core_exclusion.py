"""분석에서 제외하는 TMA core의 단일 출처(single source of truth).

두 종류의 제외:
  - qc_fail            : TMA core QC 미통과 (TMA1 D4, D8)
  - analysis_outlier   : QC는 통과했으나 분석 왜곡을 일으키는 outlier
                         (TMA1 D6: luminal_mature가 이 단일 코어에 82% 쏠려
                          DCIS vs MIBC 비교에서 사실상 단일 환자 효과)

tumor pool(CL3+4+10) 분석에서는 추가로 subcluster 7(Luminal A)도 제외:
  TMA1_D6 단일코어에 81.7% 집중 → effective MIBC core 1.4 (단일환자 효과).

새 분석 스크립트는 이 모듈의 keep_mask(코어) / keep_mask_str(CSV core 문자열) /
keep_mask_tumor(코어+Luminal A) 를 사용한다.
"""

import pandas as pd

EXCLUDED_CORES = {
    ("TMA1", "D4"): "qc_fail",
    ("TMA1", "D8"): "qc_fail",
    ("TMA1", "D6"): "analysis_outlier (luminal_mature 82% 단일코어 쏠림; CL3+4+10 sub7)",
}

# 결합형 core 문자열 ("SLIDE_COREID") — CSV 테이블(core 컬럼)에서 사용
EXCLUDED_CORE_STRINGS = {f"{s}_{c}" for s, c in EXCLUDED_CORES}

# tumor pool(CL3+4+10) 재분석에서 제외하는 subcluster (leiden_sub 문자열)
EXCLUDED_TUMOR_SUBCLUSTERS = {
    "7": "Luminal A (ER+/PR+): TMA1_D6 단일코어 81.7% 집중, effective MIBC core 1.4",
    "9": "Low-content: transcripts 84·genes 33·nucleus 16µm² 최저 (저품질/부분세포; Fig03b_1)",
}


def keep_mask(obs, slide_col="slide", core_col="core_id"):
    """제외 코어가 아닌 세포의 boolean mask 반환 (slide/core_id 분리형)."""
    mask = pd.Series(True, index=obs.index)
    for s, c in EXCLUDED_CORES:
        mask &= ~((obs[slide_col] == s) & (obs[core_col] == c))
    return mask


def keep_mask_str(df, core_col="core"):
    """제외 코어가 아닌 행 mask (core='TMA1_D6' 결합형 문자열)."""
    return ~df[core_col].astype(str).isin(EXCLUDED_CORE_STRINGS)


def keep_mask_tumor(obs, slide_col="slide", core_col="core_id", sub_col="leiden_sub"):
    """제외 코어 + 제외 subcluster(Luminal A) 아닌 세포 mask (tumor pool 분석용)."""
    mask = keep_mask(obs, slide_col, core_col)
    if sub_col in obs.columns:
        mask &= ~obs[sub_col].astype(str).isin(EXCLUDED_TUMOR_SUBCLUSTERS)
    return mask


def filter_obs(obs, slide_col="slide", core_col="core_id"):
    """제외 코어 세포를 뺀 obs 반환."""
    return obs[keep_mask(obs, slide_col, core_col)]
