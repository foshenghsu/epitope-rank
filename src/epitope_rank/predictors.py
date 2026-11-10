"""Binding predictors.

Two implementations share one interface:

* `MHCflurryPredictor` - the real thing. MHCflurry 2.x runs locally, installs as a
  Python package and needs no licence server, which is why it is the default here.
* `ToyPredictor` - a deterministic stand-in with no ML in it at all. It exists so the
  test suite and CI can exercise the pipeline without downloading models. Its numbers
  are NOT predictions and must never be reported as such.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol

import pandas as pd

# Columns every predictor must return, one row per (peptide, allele).
SCHEMA = ("peptide", "allele", "affinity_nm", "percentile")


class BindingPredictor(Protocol):
    """Scores peptides against HLA alleles."""

    name: str

    def predict(self, peptides: Sequence[str], alleles: Sequence[str]) -> pd.DataFrame:
        """Return a long-format frame with the columns in `SCHEMA`.

        `affinity_nm` is predicted binding affinity in nanomolar, lower being
        tighter. `percentile` is the peptide's rank against a background of random
        peptides for that allele, lower being stronger; it is the measure to
        threshold on, because raw affinity is not comparable between alleles.
        """
        ...


class MHCflurryPredictor:
    """Adapter over `mhcflurry.Class1PresentationPredictor`.

    Imported lazily so that installing this package without the `mhcflurry` extra
    still lets everything else run.
    """

    name = "mhcflurry"

    def __init__(self) -> None:
        try:
            from mhcflurry import Class1PresentationPredictor
        except ImportError as exc:  # pragma: no cover - depends on optional extra
            raise ImportError(
                "MHCflurry is not installed. Install it with:\n"
                "    pip install 'epitope-rank[mhcflurry]'\n"
                "    mhcflurry-downloads fetch models_class1_presentation"
            ) from exc

        try:
            self._predictor = Class1PresentationPredictor.load()
        except Exception as exc:  # pragma: no cover - depends on local model files
            raise RuntimeError(
                "MHCflurry is installed but its models could not be loaded. Run:\n"
                "    mhcflurry-downloads fetch models_class1_presentation"
            ) from exc

    def predict(self, peptides: Sequence[str], alleles: Sequence[str]) -> pd.DataFrame:
        if not peptides:
            return pd.DataFrame(columns=list(SCHEMA))

        raw = self._predictor.predict(
            peptides=list(peptides),
            alleles={a: [a] for a in alleles},
            verbose=0,
        )

        # MHCflurry's column names have moved between releases, so resolve them
        # rather than assuming. Fail with a readable message if the shape changed.
        allele_col = _first_present(raw, ["best_allele", "sample_name", "allele"])
        affinity_col = _first_present(raw, ["affinity", "mhcflurry_affinity"])
        pct_col = _first_present(
            raw,
            [
                "presentation_percentile",
                "affinity_percentile",
                "mhcflurry_affinity_percentile",
            ],
        )

        out = pd.DataFrame(
            {
                "peptide": raw["peptide"].astype(str),
                "allele": raw[allele_col].astype(str),
                "affinity_nm": pd.to_numeric(raw[affinity_col], errors="coerce"),
                "percentile": pd.to_numeric(raw[pct_col], errors="coerce"),
            }
        )
        return out.loc[:, list(SCHEMA)]


def _first_present(frame: pd.DataFrame, candidates: Sequence[str]) -> str:
    for col in candidates:
        if col in frame.columns:
            return col
    raise RuntimeError(
        f"MHCflurry returned none of the expected columns {list(candidates)}; "
        f"got {list(frame.columns)}. Pin a supported MHCflurry version."
    )


class ToyPredictor:
    """Deterministic placeholder scores for tests and offline demos.

    The score is a crude hydrophobicity-and-anchor heuristic. It is reproducible and
    has the right shape for the pipeline to consume, and that is its only purpose.
    """

    name = "toy"

    # Kyte-Doolittle hydropathy, used here only to make the numbers vary in a
    # stable, explainable way.
    _HYDROPATHY = {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
        "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
        "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
        "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
    }
    # Residues favoured at the C-terminal anchor by several common class I alleles.
    _C_ANCHOR = frozenset("VLIMFYK")

    def predict(self, peptides: Sequence[str], alleles: Sequence[str]) -> pd.DataFrame:
        rows = []
        n = len(peptides)
        for allele in alleles:
            # Allele-specific offset so the panel does not return identical rankings.
            offset = (sum(ord(c) for c in allele) % 7) * 0.1
            strengths = []
            for pep in peptides:
                mean_hydro = sum(self._HYDROPATHY.get(c, 0.0) for c in pep) / max(len(pep), 1)
                anchor_bonus = 1.5 if pep[-1] in self._C_ANCHOR else 0.0
                strengths.append((mean_hydro + anchor_bonus - offset, pep))

            # Rank highest strength first, then derive the reported numbers from the
            # rank rather than from the raw heuristic. A percentile IS a rank against
            # a background, so this at least has the right meaning, and it keeps the
            # output on a realistic scale for any input. Ties break on the peptide
            # sequence, which keeps the whole thing deterministic.
            strengths.sort(key=lambda t: (-t[0], t[1]))
            for i, (_, pep) in enumerate(strengths):
                fraction = (i + 1) / n  # (0, 1]
                percentile = round(min(100.0, max(0.01, fraction * 100.0)), 3)
                affinity = round(max(1.0, 20000.0 * fraction * fraction), 2)
                rows.append(
                    {
                        "peptide": pep,
                        "allele": allele,
                        "affinity_nm": affinity,
                        "percentile": percentile,
                    }
                )
        return pd.DataFrame(rows, columns=list(SCHEMA))


def get_predictor(name: str) -> BindingPredictor:
    """Build a predictor by name (`mhcflurry` or `toy`)."""
    key = name.strip().lower()
    if key == "mhcflurry":
        return MHCflurryPredictor()
    if key == "toy":
        return ToyPredictor()
    raise ValueError(f"Unknown predictor {name!r}. Choose 'mhcflurry' or 'toy'.")
