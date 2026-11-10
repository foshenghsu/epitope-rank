# epitope-rank

Tile a protein antigen into peptides, score them against an HLA class I panel, and get back a ranked shortlist of candidate epitopes — the first decision in a T-cell epitope campaign, written down so it can be re-run and argued with.

The question it answers: *given this antigen and the HLA alleles I care about, which peptides should I synthesise first?*

```bash
epitope-rank \
  --fasta data/NY-ESO-1_P78358.fasta \
  --alleles 'HLA-A*02:01,HLA-A*01:01,HLA-B*07:02,HLA-B*08:01' \
  --lengths 9,10 \
  --top 10 --plot --out results
```

## Install

```bash
pip install -e '.[plot]'
```

The real predictor is an optional extra, because it pulls in TensorFlow and needs a one-time model download:

```bash
pip install -e '.[mhcflurry,plot]'
mhcflurry-downloads fetch models_class1_presentation
```

[MHCflurry](https://github.com/openvax/mhcflurry) is the default because it runs locally as a Python package with no licence server, so a pipeline built on it stays reproducible. Without it, `--predictor toy` runs the same pipeline on a deterministic placeholder (see *Limitations*).

## How it ranks

Two things decide whether a peptide is worth making, and they are kept separate rather than mashed into one score:

- **strength** — how tightly the best allele in the panel binds it, as a percentile rank against a background of random peptides. Percentile, not raw nM, because affinity is not comparable between alleles.
- **breadth** — how many alleles in the panel bind it at all. A peptide that binds four alleles covers more of the patient population per peptide synthesised, which is the thing that actually costs money.

A peptide must clear both cutoffs (`--affinity-nm`, default 500 nM, and `--percentile`, default 2.0) to count as bound by an allele. `--sort` then decides which property leads:

- `breadth` (default) — for an unselected cohort, where each peptide synthesised should cover as much of the panel as possible.
- `strength` — for a cohort already selected for one allele, where a potent restricted epitope is the goal and promiscuous mid-affinity binders are noise.

Defaults are the conventions in the class I literature, not truths — they are flags so a campaign can tighten them on purpose.

Output lands in `--out`:

| file | contents |
|---|---|
| `shortlist.tsv` | ranked candidates, one row per peptide position |
| `all_scores.tsv` | every peptide × allele pair, with an `is_binder` column |
| `landscape.png` | best percentile against position, shortlist circled (`--plot`) |

`all_scores.tsv` is kept deliberately: when someone asks why a peptide *didn't* make the list, the answer should be a row, not a shrug.

## Worked example, and what it shows

NY-ESO-1 (UniProt [P78358](https://rest.uniprot.org/uniprotkb/P78358.fasta), 180 aa) is a good test antigen because part of the answer is already known: it has a well-described immunodominant HLA-A\*02:01 epitope at residues **157–165, SLLMWITQC**.

The command at the top of this README (MHCflurry 2.1, four alleles, 343 peptides) shortlists six:

| rank | peptide | start | best allele | nM | %ile | breadth |
|---|---|---|---|---|---|---|
| 1 | AADHRQLQL | 139 | HLA-B\*07:02 | 248.8 | 0.39 | 3/4 |
| 2 | APRGPHGGAA | 60 | HLA-B\*07:02 | 16.0 | 0.14 | 2/4 |
| 3 | GPHGGAASGL | 63 | HLA-B\*07:02 | 24.5 | 0.19 | 2/4 |
| 4 | APPLPVPGV | 113 | HLA-B\*07:02 | 75.7 | 0.47 | 2/4 |
| 5 | LAMPFATPM | 92 | HLA-B\*07:02 | 151.5 | 0.65 | 2/4 |
| 6 | EAELARRSL | 101 | HLA-B\*07:02 | 219.9 | 0.88 | 2/4 |

**The known epitope is not in that list**, and the reason is the interesting part. `SLLMWITQC` is called a binder for HLA-A\*02:01 (375 nM, 1.64th percentile) and a non-binder for the other three alleles — so breadth-first ranking, which rewards covering more of the panel per peptide synthesised, puts six promiscuous B\*07:02 binders above it.

Sort by strength instead and it appears, but barely:

```bash
epitope-rank --fasta data/NY-ESO-1_P78358.fasta \
  --alleles 'HLA-A*02:01' --lengths 9 --top 0 --sort strength
```

15 peptides clear the cutoffs; `SLLMWITQC` ranks **14th of 15**, while its overlapping neighbour `LMWITQCFL` (159–167) ranks 2nd at 21 nM.

Two things are worth taking from that, and both are the point of writing this down rather than eyeballing a predictor's web form:

1. **Ranking choice changes the answer more than the predictor does.** Breadth and strength are different campaigns, so `--sort` is a flag rather than a buried default.
2. **Predicted affinity is not immunodominance.** The canonical A\*02:01 epitope of this antigen is a mid-affinity binder that a strength sort nearly discards. A pipeline that quietly took the top 10 by affinity would have missed it, which is an argument for keeping `all_scores.tsv` and for treating the shortlist as a synthesis queue rather than an answer.

Running with `--predictor toy` reproduces the mechanics without model downloads. Its numbers are a placeholder, not predictions, and every such run says so on stderr.

## Limitations

Worth stating plainly, since a tool like this is easy to over-trust:

- **Binding is not immunogenicity.** Predicted MHC binding is one necessary condition. It says nothing about proteasomal processing, TAP transport, T-cell repertoire, tolerance, or whether a response is protective.
- **Class I only.** 8–11mers in the groove. Class II is a different problem and is not modelled.
- **The `toy` predictor is not a predictor.** It is a hydrophobicity-and-anchor heuristic whose percentile is a rank against the peptides in that run. It exists so tests and CI can exercise the pipeline without downloading models, and every run that uses it prints a warning. Do not report its numbers.
- **No self-peptide filtering.** Candidates are not checked against the human proteome for cross-reactivity, which real campaigns must do.
- **No allele frequency weighting.** Breadth counts alleles in your panel equally; population coverage depends on frequencies in the population you are treating.
- **Non-canonical residues are rejected**, not scored around. A sequence containing `X` fails loudly, because predictors return meaningless numbers for unknown residues.

## Development

```bash
pip install -e '.[dev,plot]'
pytest          # 23 tests
ruff check .
```

Tests use the deterministic predictor and a stub, so the suite needs no model downloads and no network. CI runs the same thing on Python 3.10–3.13.

## Licence

MIT. The example sequence is from UniProt (CC BY 4.0).
