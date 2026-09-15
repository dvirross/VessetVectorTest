# VessetVectorTest

**Null-Model Choice for Rule-Defined Pattern Counts in Clustered Longitudinal Sequences:
A Case Study of Halachic Menstrual Anticipation Rules**

Dvir Ross — Department of Software Engineering, Shenkar College of Engineering, Design and Art, Ramat Gan, Israel;
Department of Computer Science, SCE — Shamoon College of Engineering, Beer Sheva, Israel

Target: *Biometrical Journal* (Case Study).

---

## Overview

Halachic *vesset* rules define exact arithmetic events in a woman's sequence of
menstrual intervals (three equal consecutive intervals, a constant first or second
difference, a repeated weekly anchor). This repository asks whether such events occur
in 1,554 cycles from 118 women (Fehring et al., 2013) more often than chance, and shows
that the answer depends on which exchangeability null defines "chance":

| Null | What is randomised | Joint Mahalanobis test |
|---|---|---|
| `H_G` global permutation | all cycle lengths across women | D_M = 5.08, p < .001 |
| `H_iid` multinomial | i.i.d. draws from the pooled distribution | D_M = 4.87, p < .001 |
| `H_W` within-woman permutation | order of each woman's own cycles only | D_M = 2.93, p = .13 |

Under `H_W` no pattern is individually significant after Holm correction (smallest
Holm-adjusted p = .19) and all observed counts fall at or below their null means. The
global signal is driven by Haflaga (26 observed vs 11.3 expected under `H_G`, z = 4.5).

A simulation study checks the empirical size of every test (including a split-batch
check of the covariance estimate used by the joint statistic), the power of the `H_W`
tests against AR(1) and exact-repetition alternatives, and the sensitivity of the
conclusions to unequal follow-up.

## Repository structure

```
VessetVectorTest/
├── data/FilteredData.csv        1,554 cycles, 118 women; CycleNumber contiguous within every woman
├── scripts/
│   ├── patterns.py              counting rules (vectorised + reference loop; equivalence test), null generators,
│   │                            two-sided Monte Carlo p, Holm, Mahalanobis test
│   ├── rerun_nulls.py           B = 50,000 replicates under H_G, H_iid, H_W (seed 17) -> results/null_replicates.npz
│   ├── sim_validation.py        size calibration, power vs alternatives, truncation -> results/sim_validation.npz
│   ├── postprocess.py           every number reported in the paper -> results/summary.json
│   └── build_notebook.py        generates notebooks/analysis.ipynb
├── results/                     cached Monte Carlo output (see above)
├── notebooks/analysis.ipynb     single source of truth for all analyses and figures (fig1, fig3-fig10, figA1)
├── paper/
│   ├── vesset_stat.tex          manuscript (LaTeX, author-year references)
│   ├── references.bib
│   ├── cover_letter.tex
│   └── figures/                 PDF + PNG figures; gen_consort_flow.py generates fig2
├── revision_reports/            reviewer reports, author responses, triage matrices
├── sim_results.npz, strat_full_118.pkl   legacy replicate files (superseded; computed on 1,553 cycles)
└── requirements.txt
```

## Data

The dataset was collected by Fehring et al. (2013) in a randomised trial comparing two
internet-supported natural family planning methods and is archived at Marquette
University: https://epublications.marquette.edu/data_nfp/7/

`data/FilteredData.csv` is the filtered file used in the analysis (inclusion criteria:
cycle length 18–54 days inclusive; complete cycle-length record). SHA-256:
`508ee88efb18bcd29c7ed6f841827377bb8d72a8c7dd38fdeb40cf08691f7811`. Every retained
woman's `CycleNumber` runs 1, 2, …, n_i without gaps, so no pattern window bridges an
excluded cycle; `scripts/patterns.load_data` asserts this.

**Reproducibility boundary.** All analyses in the paper are reproducible from this file
and the scripts below. The original script that derived the filtered file from the raw
Marquette archive was not preserved, so that preprocessing step cannot be reproduced
exactly, and the raw women/cycle counts before filtering are not available here.

## Reproducing the analysis

```bash
pip install -r requirements.txt
python scripts/patterns.py data/FilteredData.csv      # verifies counting code, prints observed counts
python scripts/rerun_nulls.py 50000 17                # ~1 min: replicates under the three nulls
python scripts/sim_validation.py                      # ~20 min on 4 cores: size / power / truncation
python scripts/postprocess.py                         # all reported statistics -> results/summary.json
python scripts/build_notebook.py
(cd notebooks && jupyter nbconvert --to notebook --execute --inplace analysis.ipynb)   # regenerates figures
python paper/figures/gen_consort_flow.py              # fig2
```

Pattern order everywhere: Haflaga, Dilug, Week, Week-Dilug, Dilug-in-Dilug.
Observed counts: 26, 61, 32, 18, 27.

### Compiling the paper

```bash
cd paper
pdflatex -interaction=nonstopmode vesset_stat.tex
bibtex vesset_stat
pdflatex -interaction=nonstopmode vesset_stat.tex
pdflatex -interaction=nonstopmode vesset_stat.tex
```

Known non-fatal warnings: `titlesec` "entered in horizontal mode" (from the run-in
`\paragraph` format).

## Citation

If you use this code or data, please cite the manuscript above and the original dataset
(Fehring, R. J., Schneider, M., & Barron, K. (2013). Randomized comparison of two
Internet-supported methods of natural family planning. *Contraception*, 88(1), 24–30).
