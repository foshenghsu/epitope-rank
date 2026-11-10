"""Turn per-allele binding scores into one ranked shortlist.

The ranking answers a practical question: given a panel of HLA alleles that covers
the patient population you care about, which peptides are worth synthesising first?
Two properties matter and are kept separate in the output:

* strength - how tightly the best allele in the panel binds the peptide;
* breadth  - how many alleles in the panel bind it at all. A peptide that binds
  several alleles covers more of the population per peptide synthesised.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from epitope_rank.predictors import BindingPredictor
from epitope_rank.tiling import Peptide, unique_sequences

# 500 nM is the convention for "a binder" in the class I literature, and 2.0 is the
# usual strong-binder percentile cut. Both are defaults, not truths - they are
# exposed so a campaign can tighten or loosen them deliberately.
DEFAULT_AFFINITY_NM = 500.0
DEFAULT_PERCENTILE = 2.0


SORT_CHOICES = ("breadth", "strength")


@dataclass(frozen=True)
class RankingConfig:
    """Thresholds, sort order and output size for one ranking run."""

    affinity_nm: float = DEFAULT_AFFINITY_NM
    percentile: float = DEFAULT_PERCENTILE
    top_n: int | None = 25
    sort_by: str = "breadth"

    def __post_init__(self) -> None:
        if self.affinity_nm <= 0:
            raise ValueError("affinity_nm must be positive")
        if not 0 < self.percentile <= 100:
            raise ValueError("percentile must be in (0, 100]")
        if self.top_n is not None and self.top_n < 1:
            raise ValueError("top_n must be at least 1, or None for no limit")
        if self.sort_by not in SORT_CHOICES:
            raise ValueError(f"sort_by must be one of {SORT_CHOICES}, got {self.sort_by!r}")


def score_peptides(
    peptides: Sequence[Peptide],
    alleles: Sequence[str],
    predictor: BindingPredictor,
) -> pd.DataFrame:
    """Score every distinct peptide against every allele.

    Returns the long-format predictor output joined to peptide positions, so a
    peptide occurring twice in the antigen appears once per position.
    """
    if not alleles:
        raise ValueError("No HLA alleles given")
    if not peptides:
        return pd.DataFrame(columns=["peptide", "allele", "affinity_nm", "percentile",
                                     "start", "end", "length", "protein"])

    scores = predictor.predict(unique_sequences(peptides), list(alleles))
    scores = scores.assign(
        affinity_nm=scores["affinity_nm"].round(1),
        percentile=scores["percentile"].round(3),
    )
    positions = pd.DataFrame(
        {
            "peptide": [p.sequence for p in peptides],
            "start": [p.start for p in peptides],
            "end": [p.end for p in peptides],
            "length": [p.length for p in peptides],
            "protein": [p.protein for p in peptides],
        }
    ).drop_duplicates(subset=["peptide", "start"])

    return scores.merge(positions, on="peptide", how="inner")


def rank_peptides(
    peptides: Sequence[Peptide],
    alleles: Sequence[str],
    predictor: BindingPredictor,
    config: RankingConfig | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score, then collapse to one row per peptide position, ranked.

    Returns `(shortlist, all_scores)`. `all_scores` is every peptide-allele pair, kept
    so a rejected peptide can be explained rather than just missing.

    Shortlist columns:
        peptide, start, end, length, protein,
        best_allele, best_affinity_nm, best_percentile,
        n_alleles_bound, alleles_bound, rank
    """
    cfg = config or RankingConfig()
    all_scores = score_peptides(peptides, alleles, predictor)
    if all_scores.empty:
        return all_scores.copy(), all_scores

    binder = (all_scores["affinity_nm"] <= cfg.affinity_nm) & (
        all_scores["percentile"] <= cfg.percentile
    )
    all_scores = all_scores.assign(is_binder=binder)

    group_cols = ["peptide", "start", "end", "length", "protein"]
    # The reported "best" allele must be one that actually binds. Because the
    # cutoff is an AND of affinity and percentile, the lowest-percentile allele
    # can fail the nM cut while a different allele passes both - reporting the
    # former would contradict `alleles_bound` in the same row. Sorting binders
    # first keeps the two columns consistent. drop_duplicates (rather than
    # groupby().first(), which takes the first non-null value of each column
    # independently) keeps allele, affinity and percentile from the same row.
    strongest = (
        all_scores.sort_values(["is_binder", "percentile"], ascending=[False, True])
        .drop_duplicates(subset=group_cols, keep="first")
        .reset_index(drop=True)
    )
    strongest = strongest.rename(
        columns={
            "allele": "best_allele",
            "affinity_nm": "best_affinity_nm",
            "percentile": "best_percentile",
        }
    )

    bound = all_scores[all_scores["is_binder"]]
    breadth = (
        bound.groupby(group_cols, as_index=False)
        .agg(
            n_alleles_bound=("allele", "nunique"),
            alleles_bound=("allele", lambda s: ";".join(sorted(set(s)))),
        )
    )

    shortlist = strongest.merge(breadth, on=group_cols, how="left")
    shortlist["n_alleles_bound"] = shortlist["n_alleles_bound"].fillna(0).astype(int)
    shortlist["alleles_bound"] = shortlist["alleles_bound"].fillna("")

    # Keep only peptides at least one allele binds.
    shortlist = shortlist[shortlist["n_alleles_bound"] > 0]

    # Which of breadth and strength leads the sort is a real decision, not a
    # default to bury. Breadth suits an unselected cohort, where covering more
    # alleles per peptide synthesised is the win. Strength suits a cohort already
    # selected for one allele, where a potent allele-restricted epitope is the
    # point and promiscuous mid-affinity binders are noise.
    if cfg.sort_by == "breadth":
        by, ascending = ["n_alleles_bound", "best_percentile", "start"], [False, True, True]
    else:
        by, ascending = ["best_percentile", "n_alleles_bound", "start"], [True, False, True]
    shortlist = shortlist.sort_values(by, ascending=ascending).reset_index(drop=True)
    shortlist.insert(0, "rank", shortlist.index + 1)

    # Round for readability. Full float precision in a shortlist implies a
    # confidence these predictions do not have.
    shortlist["best_affinity_nm"] = shortlist["best_affinity_nm"].round(1)
    shortlist["best_percentile"] = shortlist["best_percentile"].round(3)

    if cfg.top_n is not None:
        shortlist = shortlist.head(cfg.top_n).copy()

    columns = [
        "rank", "peptide", "start", "end", "length", "protein",
        "best_allele", "best_affinity_nm", "best_percentile",
        "n_alleles_bound", "alleles_bound",
    ]
    return shortlist.loc[:, columns], all_scores
