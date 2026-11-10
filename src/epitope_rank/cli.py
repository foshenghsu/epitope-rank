"""Command-line entry point: FASTA in, ranked shortlist out."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from epitope_rank import __version__
from epitope_rank.fasta import display_name, read_fasta
from epitope_rank.predictors import get_predictor
from epitope_rank.ranking import (
    DEFAULT_AFFINITY_NM,
    DEFAULT_PERCENTILE,
    SORT_CHOICES,
    RankingConfig,
    rank_peptides,
)
from epitope_rank.tiling import DEFAULT_LENGTHS, tile


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="epitope-rank",
        description=(
            "Tile a protein antigen into peptides, score them against an HLA class I "
            "panel, and write a ranked shortlist of candidate epitopes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fasta", required=True, type=Path, help="Antigen FASTA file")
    p.add_argument(
        "--alleles",
        required=True,
        help="Comma-separated HLA alleles, e.g. 'HLA-A*02:01,HLA-B*07:02'",
    )
    p.add_argument(
        "--lengths",
        default=",".join(str(n) for n in DEFAULT_LENGTHS),
        help="Comma-separated peptide lengths to tile",
    )
    p.add_argument(
        "--predictor",
        default="mhcflurry",
        choices=["mhcflurry", "toy"],
        help="'toy' is a deterministic stand-in for offline runs, not a prediction",
    )
    p.add_argument("--affinity-nm", type=float, default=DEFAULT_AFFINITY_NM,
                   help="Binder cutoff in nM (lower is stricter)")
    p.add_argument("--percentile", type=float, default=DEFAULT_PERCENTILE,
                   help="Binder cutoff on percentile rank (lower is stricter)")
    p.add_argument(
        "--sort",
        default="breadth",
        choices=list(SORT_CHOICES),
        help=(
            "'breadth' ranks peptides bound by more alleles first (unselected cohort); "
            "'strength' ranks the tightest binder first (single-allele cohort)"
        ),
    )
    p.add_argument("--top", type=int, default=25, help="Rows in the shortlist; 0 for all")
    p.add_argument(
        "--name",
        default=None,
        help="Antigen name for output rows and the figure title "
        "(default: derived from the FASTA header)",
    )
    p.add_argument("--out", type=Path, default=Path("results"), help="Output directory")
    p.add_argument("--plot", action="store_true", help="Also write landscape.png")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def _split(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    records = read_fasta(args.fasta)
    if len(records) > 1:
        print(
            f"note: {len(records)} records in {args.fasta.name}; scoring the first "
            f"({records[0].name})",
            file=sys.stderr,
        )
    record = records[0]

    alleles = _split(args.alleles)
    lengths = [int(n) for n in _split(args.lengths)]
    antigen = args.name or display_name(record.name)
    peptides = tile(record.sequence, lengths=lengths, protein=antigen)
    if not peptides:
        print(
            f"error: antigen is {len(record)} aa, shorter than every requested "
            f"peptide length {sorted(lengths)}",
            file=sys.stderr,
        )
        return 2

    predictor = get_predictor(args.predictor)
    config = RankingConfig(
        affinity_nm=args.affinity_nm,
        percentile=args.percentile,
        top_n=None if args.top == 0 else args.top,
        sort_by=args.sort,
    )
    shortlist, all_scores = rank_peptides(peptides, alleles, predictor, config)

    args.out.mkdir(parents=True, exist_ok=True)
    shortlist_path = args.out / "shortlist.tsv"
    scores_path = args.out / "all_scores.tsv"
    shortlist.to_csv(shortlist_path, sep="\t", index=False)
    all_scores.to_csv(scores_path, sep="\t", index=False)

    print(
        f"{record.name}: {len(record)} aa -> {len(peptides)} peptides "
        f"x {len(alleles)} alleles ({predictor.name})"
    )
    print(
        f"{len(shortlist)} shortlisted at <= {config.affinity_nm:g} nM and "
        f"<= {config.percentile:g} percentile, sorted by {config.sort_by}"
    )
    print(f"wrote {shortlist_path}")
    print(f"wrote {scores_path}")

    if args.plot:
        from epitope_rank.plotting import plot_landscape

        png = plot_landscape(
            all_scores, shortlist, args.out / "landscape.png",
            protein_length=len(record), title=f"{antigen} ({predictor.name})",
        )
        print(f"wrote {png}")

    if args.predictor == "toy":
        print(
            "\nWARNING: the 'toy' predictor is a deterministic placeholder. Its numbers "
            "are not binding predictions. Use --predictor mhcflurry for real scores.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
