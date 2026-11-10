import pytest

from epitope_rank.fasta import parse_fasta


def test_parses_multiline_and_upper_cases():
    records = parse_fasta(">test antigen\nacdef\nghikl\n")
    assert len(records) == 1
    assert records[0].name == "test antigen"
    assert records[0].sequence == "ACDEFGHIKL"
    assert len(records[0]) == 10


def test_parses_multiple_records():
    assert len(parse_fasta(">a\nACDEF\n>b\nGHIKL\n")) == 2


def test_rejects_non_canonical_residues():
    # X would otherwise be handed to a predictor, which returns a meaningless score.
    with pytest.raises(ValueError, match="non-canonical"):
        parse_fasta(">bad\nACDEFX\n")


def test_rejects_empty_sequence():
    with pytest.raises(ValueError, match="no sequence"):
        parse_fasta(">empty\n")


def test_rejects_sequence_before_header():
    with pytest.raises(ValueError, match="before any"):
        parse_fasta("ACDEF\n")


def test_rejects_empty_input():
    with pytest.raises(ValueError, match="No FASTA records"):
        parse_fasta("\n\n")
