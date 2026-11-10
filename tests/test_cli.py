import pandas as pd
import pytest

from epitope_rank.cli import main

FIXTURE = ">TESTANTIGEN synthetic\nACDEFGHIKLMNPQRSTVWYACDEFGHIKL\n"


def _write(tmp_path, text=FIXTURE):
    path = tmp_path / "antigen.fasta"
    path.write_text(text, encoding="utf-8")
    return path


def test_cli_writes_shortlist_and_scores(tmp_path, capsys):
    fasta = _write(tmp_path)
    out = tmp_path / "results"
    code = main(
        [
            "--fasta", str(fasta),
            "--alleles", "HLA-A*02:01,HLA-B*07:02",
            "--predictor", "toy",
            "--percentile", "100",
            "--affinity-nm", "1e9",
            "--out", str(out),
        ]
    )
    assert code == 0

    shortlist = pd.read_csv(out / "shortlist.tsv", sep="\t")
    assert not shortlist.empty
    assert list(shortlist["rank"]) == list(range(1, len(shortlist) + 1))
    assert (out / "all_scores.tsv").exists()

    captured = capsys.readouterr()
    assert "TESTANTIGEN" in captured.out
    # The placeholder predictor must announce itself on every run.
    assert "not binding predictions" in captured.err


def test_cli_errors_when_antigen_shorter_than_peptides(tmp_path):
    fasta = _write(tmp_path, ">short\nACDEF\n")
    code = main(
        ["--fasta", str(fasta), "--alleles", "HLA-A*02:01", "--predictor", "toy",
         "--out", str(tmp_path / "r")]
    )
    assert code == 2


def test_cli_rejects_bad_residues(tmp_path):
    fasta = _write(tmp_path, ">bad\nACDEFGHIKXXX\n")
    with pytest.raises(ValueError, match="non-canonical"):
        main(["--fasta", str(fasta), "--alleles", "HLA-A*02:01", "--predictor", "toy",
              "--out", str(tmp_path / "r")])


def test_cli_plot_is_written_when_requested(tmp_path):
    pytest.importorskip("matplotlib")
    fasta = _write(tmp_path)
    out = tmp_path / "results"
    code = main(
        ["--fasta", str(fasta), "--alleles", "HLA-A*02:01", "--predictor", "toy",
         "--percentile", "100", "--affinity-nm", "1e9", "--out", str(out), "--plot"]
    )
    assert code == 0
    png = out / "landscape.png"
    assert png.exists() and png.stat().st_size > 0
