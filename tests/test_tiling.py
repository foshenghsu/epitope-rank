from epitope_rank.tiling import tile, unique_sequences


def test_tiling_counts_and_positions():
    seq = "ACDEFGHIKLM"  # 11 residues
    peptides = tile(seq, lengths=[9])
    assert len(peptides) == 3  # 11 - 9 + 1
    assert peptides[0].sequence == "ACDEFGHIK"
    assert peptides[0].start == 1
    assert peptides[0].end == 9
    assert peptides[-1].start == 3
    assert peptides[-1].sequence == seq[2:11]


def test_tiling_multiple_lengths_sorted_by_position():
    peptides = tile("ACDEFGHIKLM", lengths=[10, 9])
    starts = [p.start for p in peptides]
    assert starts == sorted(starts)
    assert {p.length for p in peptides} == {9, 10}


def test_length_longer_than_protein_is_skipped():
    assert tile("ACDEF", lengths=[9]) == []


def test_duplicate_lengths_do_not_duplicate_peptides():
    assert len(tile("ACDEFGHIKLM", lengths=[9, 9, 9])) == 3


def test_unique_sequences_dedupes_but_keeps_order():
    # The repeat means AAAAAAAAA is tiled at more than one start position.
    peptides = tile("AAAAAAAAAACDEF", lengths=[9])
    uniq = unique_sequences(peptides)
    assert len(uniq) == len({p.sequence for p in peptides})
    assert uniq[0] == peptides[0].sequence
