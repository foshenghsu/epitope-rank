"""Minimal FASTA reading, with enough validation to fail loudly on bad input."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# The 20 canonical amino acids. Anything else (X, B, Z, U, O, gaps) is rejected
# rather than silently scored, because predictors give meaningless numbers for
# peptides containing unknown residues.
CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class Record:
    """One FASTA entry."""

    name: str
    sequence: str

    def __len__(self) -> int:
        return len(self.sequence)


def parse_fasta(text: str) -> list[Record]:
    """Parse FASTA text into records, upper-casing and validating each sequence."""
    records: list[Record] = []
    name: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        if name is None:
            return
        seq = "".join(chunks).upper()
        if not seq:
            raise ValueError(f"FASTA record {name!r} has no sequence")
        bad = sorted(set(seq) - CANONICAL_AA)
        if bad:
            raise ValueError(
                f"FASTA record {name!r} contains non-canonical residues: {''.join(bad)}. "
                "Remove or substitute them before scoring."
            )
        records.append(Record(name=name, sequence=seq))

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            name = line[1:].strip() or "unnamed"
            chunks = []
        else:
            if name is None:
                raise ValueError("FASTA data began before any '>' header line")
            chunks.append(line.replace(" ", ""))
    flush()

    if not records:
        raise ValueError("No FASTA records found")
    return records


def read_fasta(path: str | Path) -> list[Record]:
    """Read and validate a FASTA file."""
    return parse_fasta(Path(path).read_text(encoding="utf-8"))
