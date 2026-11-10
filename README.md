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

- **strength** — how tightly the best allele in the panel binds it, as an *affinity* percentile rank against a background of random peptides. Percentile, not raw nM, because affinity is not comparable between alleles. (MHCflurry also reports a `presentation_percentile`, which folds in antigen processing; this tool asks for the affinity percentile specifically, because the 2.0 cutoff below is an affinity-percentile convention and applying it to a presentation score would be a different filter.)
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

The command at the top of this README, run with `--out examples/results` (MHCflurry 2.3.25, four alleles, 343 peptides, 41 of which bind at least one allele), shortlists ten — the numbers below are copied from the committed [`examples/results/shortlist.tsv`](examples/results/shortlist.tsv):

| rank | peptide | start | best allele | nM | %ile | breadth |
|---|---|---|---|---|---|---|
| 1 | AADHRQLQL | 139 | HLA-B\*08:01 | 383.2 | 0.614 | 3/4 |
| 2 | APRGPHGGAA | 60 | HLA-B\*07:02 | 16.0 | 0.007 | 2/4 |
| 3 | APRGPHGGA | 60 | HLA-B\*07:02 | 17.9 | 0.014 | 2/4 |
| 4 | GPHGGAASGL | 63 | HLA-B\*07:02 | 24.5 | 0.046 | 2/4 |
| 5 | APPLPVPGVL | 113 | HLA-B\*07:02 | 33.5 | 0.092 | 2/4 |
| 6 | APPLPVPGV | 113 | HLA-B\*07:02 | 75.7 | 0.301 | 2/4 |
| 7 | LAMPFATPM | 92 | HLA-B\*07:02 | 151.5 | 0.528 | 2/4 |
| 8 | EAELARRSL | 101 | HLA-B\*08:01 | 353.9 | 0.584 | 2/4 |
| 9 | GPGGGAPRGP | 55 | HLA-B\*07:02 | 306.3 | 0.777 | 2/4 |
| 10 | LLMWITQCFL | 158 | HLA-A\*02:01 | 14.6 | 0.057 | 1/4 |

**The canonical 9mer itself is not on that list**, and the reason is the interesting part. `SLLMWITQC` is a predicted binder for HLA-A\*02:01 (375.3 nM, 1.202 percentile) and a non-binder for the other three, so breadth-first ranking — which rewards covering more of the panel per peptide synthesised — puts nine peptides above it.

But the **region** is not missed. `LLMWITQCFL` (158–167) contains the epitope's core, overlaps it by eight residues, and lands at rank 10 as the single tightest binder in the table at 14.6 nM — carried there by strength despite binding only one allele. Sort by strength on an A\*02:01-only panel and the frames separate further:

```bash
epitope-rank --fasta data/NY-ESO-1_P78358.fasta \
  --alleles 'HLA-A*02:01' --lengths 9 --top 0 --sort strength
```

16 nonamers clear the cutoffs. `LMWITQCFL` (159–167) ranks **1st at 21.4 nM**; the canonical `SLLMWITQC` ranks **13th of 16**, at 375.3 nM.

Three things are worth taking from that, and they are the point of writing this down rather than eyeballing a predictor's web form:

1. **Ranking choice changes the answer more than the predictor does.** Breadth and strength are different campaigns, so `--sort` is a flag rather than a buried default.
2. **Predicted affinity is not immunodominance.** The canonical A\*02:01 epitope of this antigen is a mid-affinity binder that ranks 13th of 16 among nonamers its own allele binds. Whatever makes it immunodominant in patients, this predictor is not measuring it.
3. **The frame matters as much as the region.** A pipeline reporting only peptide identities would call rank 10 a hit and the canonical epitope a miss, when they are eight residues of the same sequence. That is an argument for keeping `all_scores.tsv` and reading positions, not just names.

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
pytest          # 24 tests
ruff check .
```

Tests use the deterministic predictor and a stub, so the suite needs no model downloads and no network. CI runs the same thing on Python 3.10–3.13.

## Licence

MIT. The example sequence is from UniProt (CC BY 4.0).
