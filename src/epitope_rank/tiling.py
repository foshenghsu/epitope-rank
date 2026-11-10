"""Tile a protein sequence into overlapping peptides of the lengths an HLA class I
groove accommodates."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

# HLA class I binds 8-11mers; 9mers dominate. Class II is longer and is out of
# scope here, which the README states.
DEFAULT_LENGTHS: tuple[int, ...] = (9, 10)


@dataclass(frozen=True)
class Peptide:
    """A peptide and where it came from.

    `start` is 1-based and inclusive, matching how positions are reported in the
    immunology literature rather than Python's 0-based slicing.
    """

    sequence: str
    start: int
    protein: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def end(self) -> int:
        return self.start + self.length - 1


def tile(
    sequence: str,
    *,
    lengths: Iterable[int] = DEFAULT_LENGTHS,
    protein: str = "antigen",
) -> list[Peptide]:
    """Return every peptide of the requested lengths, in sequence order.

    Peptides are deduplicated by (sequence, start), so repeating a length does not
    produce duplicate rows. Lengths longer than the protein are skipped rather than
    raising, which keeps short test fixtures usable.
    """
    lengths = sorted({int(n) for n in lengths})
    if not lengths:
        raise ValueError("No peptide lengths given")
    if lengths[0] < 1:
        raise ValueError(f"Peptide lengths must be positive, got {lengths[0]}")

    seq = sequence.upper()
    out: list[Peptide] = []
    for n in lengths:
        if n > len(seq):
            continue
        for i in range(len(seq) - n + 1):
            out.append(Peptide(sequence=seq[i : i + n], start=i + 1, protein=protein))

    out.sort(key=lambda p: (p.start, p.length))
    return out


def unique_sequences(peptides: Sequence[Peptide]) -> list[str]:
    """Distinct peptide sequences, order preserved.

    A peptide can occur at more than one position in a repetitive antigen, and
    predictors are expensive, so each distinct sequence is scored once and the
    result is joined back onto every position it came from.
    """
    seen: set[str] = set()
    out: list[str] = []
    for p in peptides:
        if p.sequence not in seen:
            seen.add(p.sequence)
            out.append(p.sequence)
    return out
