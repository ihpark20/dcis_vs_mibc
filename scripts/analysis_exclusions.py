"""Cores left out of the analysis, in one place.

The clustering keeps every QC-passing core, so that the cell set matches the published one.
What comes after it does not: TMA1 D6 supplies most of the cells in two of the tumor pool's
luminal-mature subclusters, and a cell state carried by a single core cannot support a
comparison between DCIS and microinvasive disease. It is therefore dropped from everything
downstream of the clustering — subclustering, boundaries, communication — and every script
that does so reads the exclusion from here.
"""

EXCLUDED_CORES = {"TMA1_D6"}
EXCLUSION_REASON = {
    "TMA1_D6": "one core supplies most of two luminal-mature tumor subclusters",
}


def drop_excluded(frame, core_column="core"):
    """Remove the excluded cores from a table carrying a core column."""
    keep = ~frame[core_column].isin(EXCLUDED_CORES)
    dropped = (~keep).sum()
    if dropped:
        names = ", ".join(sorted(set(frame.loc[~keep, core_column])))
        print(f"excluded from the analysis: {names} ({dropped:,} cells)")
    return frame[keep]
