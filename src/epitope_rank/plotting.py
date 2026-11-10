"""One figure: where the predicted binders sit along the antigen."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def plot_landscape(
    all_scores: pd.DataFrame,
    shortlist: pd.DataFrame,
    out_path: str | Path,
    *,
    protein_length: int | None = None,
    title: str | None = None,
) -> Path:
    """Draw best-percentile against start position, with shortlisted peptides marked.

    A lower percentile is a stronger binder, so the y axis is inverted: strong
    binders sit at the top, where the eye goes first. Returns the written path.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")  # headless; never try to open a window in CI
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise ImportError(
            "Plotting needs matplotlib. Install it with: pip install 'epitope-rank[plot]'"
        ) from exc

    if all_scores.empty:
        raise ValueError("Nothing to plot: no scores were produced")

    best = (
        all_scores.sort_values("percentile")
        .groupby(["peptide", "start"], as_index=False)
        .first()
    )

    fig, ax = plt.subplots(figsize=(11.5, 4.2), dpi=160)
    ax.scatter(
        best["start"], best["percentile"],
        s=14, c="#9aa3ad", alpha=0.65, linewidths=0, label="all peptides (best allele)",
    )
    if not shortlist.empty:
        ax.scatter(
            shortlist["start"], shortlist["best_percentile"],
            s=52, facecolors="none", edgecolors="#c2410c", linewidths=1.6,
            label=f"shortlist (n={len(shortlist)})", zorder=3,
        )
        # Label only the top few. A shortlist is mostly overlapping frames of the
        # same region, so neighbouring labels collide: alternating above/below is
        # not enough when two peptides share a start position. Step the offsets
        # through several distances and nudge alternate labels sideways, with a
        # leader line so a displaced label still reads as belonging to its marker.
        offsets = [(0, 14), (0, -18), (34, 24), (-34, -28), (0, 34)]
        # strict=False: a shortlist shorter than `offsets` is normal (small panels).
        for (dx, dy), row in zip(
            offsets, shortlist.head(len(offsets)).itertuples(), strict=False
        ):
            ax.annotate(
                row.peptide,
                (row.start, row.best_percentile),
                textcoords="offset points",
                xytext=(dx, dy),
                ha="center",
                va="bottom" if dy > 0 else "top",
                fontsize=7,
                color="#7c2d12",
                zorder=4,
                arrowprops=dict(
                    arrowstyle="-", color="#c2410c", linewidth=0.6,
                    shrinkA=0, shrinkB=6, alpha=0.7,
                ),
            )

    ax.set_yscale("log")
    ax.invert_yaxis()  # strong binders (low percentile) at the top
    ax.set_xlabel("peptide start position in antigen (1-based)")
    ax.set_ylabel("best-allele affinity percentile (log, lower = stronger)")
    ax.set_title(title or "Predicted binding landscape")
    if protein_length:
        ax.set_xlim(0, protein_length + 1)
    ax.grid(True, which="major", axis="both", alpha=0.18)
    # Outside the axes: in here, an unfilled legend marker reads as a data point.
    ax.legend(
        frameon=False, fontsize=8, loc="upper left",
        bbox_to_anchor=(1.01, 1.0), borderaxespad=0,
    )
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out
