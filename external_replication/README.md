# External replication — Utah Creighton Model cohort

Self-contained replication of the manuscript's analysis on an independent, openly deposited dataset
(Stanford & Najmabadi 2023, *Menstrual Cycles Length of Women in the USA and Canada, 1990–2013*,
The Hive, University of Utah, https://doi.org/10.7278/S50d-4gxs-s4hj, CC BY-NC). Nothing outside this
folder is written by the code here; counting rules are imported from `../scripts/patterns.py`.

| | |
|---|---|
| `reports/utah_replication_report.md` | all results and the side-by-side with Fehring |
| `reports/utah_preprocessing_and_deviations.md` | dataset selection, checksums, preprocessing rules, deviation log |
| `reports/utah_manuscript_subsection_draft.md` | draft subsection (not inserted in the manuscript) |
| `data/` | raw deposit (`raw/`), analysed `sequences.csv` (270 women, 2,139 cycles), `women.csv`, `preprocessing.json` |
| `results/` | B = 50,000 replicates, summary, 18–54 sensitivity, validation, comparison tables, figures, provenance |
| `scripts/` | `prepare_utah.py`, `run_external.py`, `analysis_engine.py`, `external_figures.py` |

Reproduce (from the repository root, ~30 min on 3 cores):

```bash
python external_replication/scripts/prepare_utah.py            # raw -> data/sequences.csv
python external_replication/scripts/run_external.py --nproc 3  # -> results/
```

Headline: observed counts (41, 76, 35, 40, 41); joint D_M = 9.99 / 9.53 under H_G / H_iid (p < .001),
1.20 under H_W (p = .92; LOO .86–.97). Fehring: 5.08 / 4.87 / 2.93 (p = .126).
