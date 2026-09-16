# VessetVectorTest

**Null-Model Choice for Rule-Defined Pattern Counts in Clustered Longitudinal Sequences:
A Case Study of Jewish Halachic Menstrual Anticipation Rules**

Dvir Ross — Department of Software Engineering, Shenkar College of Engineering, Design and Art, Ramat Gan, Israel;
Department of Computer Science, SCE — Shamoon College of Engineering, Beer Sheva, Israel

Target: *Biometrical Journal* (Case Study).

---

## Overview

Halachic *vesset* rules define exact arithmetic events in a woman's sequence of
menstrual intervals (three equal consecutive intervals, a constant first or second
difference, a repeated weekly anchor). This repository asks whether such events occur
in 1,554 cycles from 118 women (Fehring, 2012) more often than chance, and shows
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
├── data/RawData.csv             archived source file (Fehring, 2012): 1,665 cycles, 159 women
├── data/FilteredData.csv        1,554 cycles, 118 women; CycleNumber contiguous within every woman
├── scripts/
│   ├── filter_raw.py            data/RawData.csv -> data/FilteredData.csv (verified identical)
│   ├── Filtering_original.ipynb original filtering notebook (from the PhD repository)
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
│   ├── supplementary_for_review/  PhD dissertation (2022) and accepted B.D.D. article (2021), Hebrew; for review only
│   └── figures/                 PDF + PNG figures; gen_consort_flow.py generates fig2
├── revision_reports/            reviewer reports, author responses, triage matrices
├── sim_results.npz, strat_full_118.pkl   legacy replicate files (superseded; computed on 1,553 cycles)
└── requirements.txt
```

## Data

The dataset is the Marquette University "Menstrual Cycle Data" archive (Fehring, 2012),
collected in a randomised trial of two Internet-supported fertility-awareness-based
methods of family planning (Fehring et al., 2013):
https://epublications.marquette.edu/data_nfp/7/

`data/FilteredData.csv` is the filtered file used in the analysis. SHA-256:
`508ee88efb18bcd29c7ed6f841827377bb8d72a8c7dd38fdeb40cf08691f7811`. Every retained
woman's `CycleNumber` runs 1, 2, …, n_i without gaps, so no pattern window bridges an
excluded cycle; `scripts/patterns.load_data` asserts this.

**Preprocessing.** `scripts/filter_raw.py` derives the filtered file from the archived
source file (1,665 cycles, 159 women; all cycle lengths already 18–54 days, none
missing, so no cycle-length criterion is applied):

1. remove 3 women whose cycle records appear twice (nfp8106, nfp8109, nfp8114; 16 rows) → 1,649 cycles, 156 women
2. keep women with ≥ 5 recorded cycles (38 women, 87 cycles removed) → 1,562 cycles, 118 women
3. nfp8107's eight cycles appear as two near-identical copies (cycle 2: 27 vs 26 days); the first copy is dropped → 1,554 cycles, 118 women

```bash
python scripts/filter_raw.py data/RawData.csv     # verifies output == data/FilteredData.csv
```

`data/RawData.csv` is the author's copy of the archived source file (SHA-256
`9fa76583129db58af44adb21a0612217b207120fb61d478fb343a61291732a89`);
`scripts/Filtering_original.ipynb` is the original filtering notebook from
https://github.com/dvirross/PhD, which `filter_raw.py` re-implements.

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
(Fehring, R. J. (2012). Menstrual Cycle Data. Marquette University e-Publications,
https://epublications.marquette.edu/data_nfp/7/) and the trial article (Fehring, R. J.,
Schneider, M., Raviele, K., Rodriguez, D., & Pruszynski, J. (2013). Randomized comparison
of two Internet-supported fertility-awareness-based methods of family planning.
*Contraception*, 88(1), 24–30. https://doi.org/10.1016/j.contraception.2012.10.010).
