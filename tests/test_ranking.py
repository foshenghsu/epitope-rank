import pandas as pd
import pytest

from epitope_rank.predictors import SCHEMA, ToyPredictor
from epitope_rank.ranking import RankingConfig, rank_peptides, score_peptides
from epitope_rank.tiling import tile

ALLELES = ["HLA-A*02:01", "HLA-B*07:02"]


class StubPredictor:
    """Fixed scores, so the ranking logic is tested on its own."""

    name = "stub"

    def __init__(self, table):
        self._table = table

    def predict(self, peptides, alleles):
        rows = [r for r in self._table if r["peptide"] in set(peptides)]
        return pd.DataFrame(rows, columns=list(SCHEMA))


def test_toy_predictor_matches_schema_and_is_deterministic():
    peptides = tile("ACDEFGHIKLMNPQRST", lengths=[9])
    pred = ToyPredictor()
    first = pred.predict([p.sequence for p in peptides], ALLELES)
    second = pred.predict([p.sequence for p in peptides], ALLELES)
    assert list(first.columns) == list(SCHEMA)
    assert len(first) == len(peptides) * len(ALLELES)
    pd.testing.assert_frame_equal(first, second)


def test_scores_carry_positions():
    peptides = tile("ACDEFGHIKLMNP", lengths=[9])
    scored = score_peptides(peptides, ALLELES, ToyPredictor())
    assert {"start", "end", "length", "protein"} <= set(scored.columns)
    assert scored["start"].min() == 1


def test_ranking_prefers_breadth_then_strength():
    # BROADPEPT binds both alleles; TIGHTPEPT binds one allele harder.
    table = [
        {"peptide": "BROADPEPT", "allele": ALLELES[0], "affinity_nm": 90.0, "percentile": 0.9},
        {"peptide": "BROADPEPT", "allele": ALLELES[1], "affinity_nm": 120.0, "percentile": 1.1},
        {"peptide": "TIGHTPEPT", "allele": ALLELES[0], "affinity_nm": 5.0, "percentile": 0.05},
        {"peptide": "TIGHTPEPT", "allele": ALLELES[1], "affinity_nm": 9000.0, "percentile": 45.0},
        {"peptide": "WEAKPEPTX", "allele": ALLELES[0], "affinity_nm": 8000.0, "percentile": 30.0},
        {"peptide": "WEAKPEPTX", "allele": ALLELES[1], "affinity_nm": 9500.0, "percentile": 60.0},
    ]
    seq = "BROADPEPTTIGHTPEPTWEAKPEPTX"
    shortlist, all_scores = rank_peptides(
        tile(seq, lengths=[9]), ALLELES, StubPredictor(table)
    )

    assert list(shortlist["peptide"]) == ["BROADPEPT", "TIGHTPEPT"]
    assert shortlist.loc[0, "n_alleles_bound"] == 2
    assert shortlist.loc[0, "alleles_bound"] == ";".join(sorted(ALLELES))
    assert shortlist.loc[1, "best_allele"] == ALLELES[0]
    assert shortlist.loc[1, "best_percentile"] == 0.05
    assert list(shortlist["rank"]) == [1, 2]

    # Non-binders are explained, not silently dropped.
    assert "WEAKPEPTX" in set(all_scores["peptide"])
    assert not all_scores.loc[all_scores["peptide"] == "WEAKPEPTX", "is_binder"].any()


def test_sort_by_strength_promotes_the_tightest_binder():
    """The documented tradeoff: an allele-restricted potent binder outranks a
    promiscuous mid-affinity one only when sorting by strength."""
    table = [
        {"peptide": "BROADPEPT", "allele": ALLELES[0], "affinity_nm": 90.0, "percentile": 0.9},
        {"peptide": "BROADPEPT", "allele": ALLELES[1], "affinity_nm": 120.0, "percentile": 1.1},
        {"peptide": "TIGHTPEPT", "allele": ALLELES[0], "affinity_nm": 5.0, "percentile": 0.05},
        {"peptide": "TIGHTPEPT", "allele": ALLELES[1], "affinity_nm": 9000.0, "percentile": 45.0},
    ]
    peptides = tile("BROADPEPTTIGHTPEPT", lengths=[9])
    stub = StubPredictor(table)

    by_breadth, _ = rank_peptides(peptides, ALLELES, stub, RankingConfig(sort_by="breadth"))
    by_strength, _ = rank_peptides(peptides, ALLELES, stub, RankingConfig(sort_by="strength"))

    assert list(by_breadth["peptide"]) == ["BROADPEPT", "TIGHTPEPT"]
    assert list(by_strength["peptide"]) == ["TIGHTPEPT", "BROADPEPT"]
    assert list(by_strength["rank"]) == [1, 2]

    with pytest.raises(ValueError, match="sort_by"):
        RankingConfig(sort_by="nonsense")


def test_outputs_are_rounded_for_readability():
    table = [
        {"peptide": "LONGFLOAT", "allele": ALLELES[0],
         "affinity_nm": 248.78586484890488, "percentile": 0.3903549943855079},
    ]
    shortlist, all_scores = rank_peptides(
        tile("LONGFLOAT", lengths=[9]), ALLELES[:1], StubPredictor(table)
    )
    assert shortlist.loc[0, "best_affinity_nm"] == 248.8
    assert shortlist.loc[0, "best_percentile"] == 0.39
    assert all_scores.loc[0, "affinity_nm"] == 248.8


def test_thresholds_are_enforced():
    table = [
        {"peptide": "BORDERPEP", "allele": ALLELES[0], "affinity_nm": 400.0, "percentile": 1.5},
    ]
    peptides = tile("BORDERPEP", lengths=[9])
    stub = StubPredictor(table)

    loose, _ = rank_peptides(peptides, ALLELES[:1], stub)
    assert len(loose) == 1

    # Strict on affinity alone is enough to reject it.
    strict, _ = rank_peptides(
        peptides, ALLELES[:1], stub, RankingConfig(affinity_nm=100.0)
    )
    assert strict.empty

    # And strict on percentile alone, independently.
    strict2, _ = rank_peptides(
        peptides, ALLELES[:1], stub, RankingConfig(percentile=0.5)
    )
    assert strict2.empty


def test_top_n_limits_rows():
    peptides = tile("ACDEFGHIKLMNPQRSTVWYACDEF", lengths=[9])
    shortlist, _ = rank_peptides(
        peptides, ALLELES, ToyPredictor(), RankingConfig(percentile=100, affinity_nm=1e9, top_n=3)
    )
    assert len(shortlist) == 3


def test_empty_peptides_and_bad_config():
    empty, _ = rank_peptides([], ALLELES, ToyPredictor())
    assert empty.empty
    with pytest.raises(ValueError):
        rank_peptides(tile("ACDEFGHIK", lengths=[9]), [], ToyPredictor())
    with pytest.raises(ValueError):
        RankingConfig(percentile=0)
    with pytest.raises(ValueError):
        RankingConfig(affinity_nm=-1)
    with pytest.raises(ValueError):
        RankingConfig(top_n=0)
