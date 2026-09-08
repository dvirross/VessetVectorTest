# VessetVectorTest

**Statistical Evidence for Halachic Vesset Patterns in Natural Menstrual Cycle Data:
A Permutation and Multinomial Sampling Analysis**

Dvir Ross — Shenkar College of Engineering, Design and Art & SCE — Shamoon College of Engineering

---

## Overview

This repository contains all code, data, and paper sources for the above article.
The study tests whether halachic *vesset* patterns (rules for anticipating menstruation
in Jewish law) occur in real menstrual cycle data at rates significantly exceeding
chance, using permutation testing and a multinomial sampling null model.

**Key findings:**
- *Haflaga* (fixed interval): z = 4.51, p < .001 — only 2 of 50,000 simulations
  reached or exceeded the observed count (the simulation ceiling)
- *Dilug* (arithmetic progression): z = 1.83, p = .043
- Joint multivariate test (all 5 patterns): p < .01 under both null models
- Focused joint test (mathematically dependent sub-vector): p < .01

## Repository Structure

```
VessetVectorTest/
├── README.md
├── requirements.txt
├── sim_results.npz          # Pre-computed simulation results (50,000 × 2 iterations)
│
├── paper/
│   ├── vesset_stat.tex      # LaTeX source (compile on Overleaf)
│   ├── references.bib       # BibTeX bibliography
│   ├── vesset_stat.pdf      # Compiled PDF
│   └── figures/             # All publication figures (PDF)
│       ├── fig0_cycle_distribution.pdf
│       ├── fig1_null_distributions.pdf
│       ├── fig2_zscores.pdf
│       ├── fig3_joint_test.pdf
│       ├── fig4_cdfs.pdf
│       ├── fig5_chisq_heatmap.pdf
│       └── fig6_3d_focused.pdf
│
├── notebooks/
│   ├── analysis.ipynb          # Main analysis: permutation & multinomial tests,
│   │                           # chi-square dependence, all figures
│   └── visualization_3d.ipynb  # Interactive & static 3D joint test visualization
│
└── data/
    └── FilteredData.csv     # Filtered NFP dataset (1,554 cycles, 118 women)
```

## Data

The dataset was originally collected by Fehring et al. (2013) in a randomised trial
comparing two internet-supported natural family planning methods, and is publicly
archived at Marquette University:
https://epublications.marquette.edu/data_nfp/7/

`FilteredData.csv` is the filtered version (cycles outside 18–54 days removed;
duplicate records resolved) used in the analysis.

## Reproducing the Analysis

### Setup

```bash
git clone https://github.com/dvirross/VessetVectorTest.git
cd VessetVectorTest
pip install -r requirements.txt
```

### Running the notebooks

```bash
jupyter notebook notebooks/analysis.ipynb
```

The main notebook (`analysis.ipynb`) loads pre-computed simulation results from
`sim_results.npz` for speed. To re-run the full 50,000-iteration simulations
from scratch, replace the "Load pre-computed results" cell with the commented-out
simulation loop (takes approximately 15 minutes).

### Compiling the paper

Upload the contents of `paper/` to [Overleaf](https://overleaf.com) as a new project,
set `vesset_stat.tex` as the main file, and compile. All figures are included in
`paper/figures/`.

## Pre-computed Simulation Results

`sim_results.npz` contains the results of 50,000 permutation iterations and 50,000
multinomial sampling iterations (random seed 17), stored as NumPy arrays:

```python
import numpy as np
d = np.load('sim_results.npz')
# d['perm']     — shape (50000, 5): permutation null counts
# d['multi']    — shape (50000, 5): multinomial null counts
# d['observed'] — shape (5,): observed pattern counts [26, 61, 32, 18, 27]
```

Pattern order: Haflaga, Dilug, Week, Week-Dilug, Dilug-in-Dilug.

## Citation

If you use this code or data in your own work, please cite:

```bibtex
@article{ross2024vessetstat,
  author  = {Dvir Ross},
  title   = {Statistical Evidence for Halachic Vesset Patterns in Natural
             Menstrual Cycle Data: A Permutation and Multinomial Sampling Analysis},
  journal = {[Journal name]},
  year    = {2024},
  note    = {Under review}
}
```

The PhD dissertation on which this work is based:

```bibtex
@phdthesis{ross2022phd,
  author  = {Dvir Ross},
  title   = {Probabilistic and Statistical Analysis of the Menstrual Cycle
             in a Halachic Context},
  school  = {Ariel University in Samaria, Ariel, Israel},
  year    = {2022},
  note    = {In Hebrew. Supervisors: Prof. E. Merzbach, Dr. E. Shmerling,
             Prof. A. Domoshnitsky}
}
```

## License

Code: MIT License.
Data: subject to the terms of the original Marquette University data archive
(https://epublications.marquette.edu/data_nfp/7/).

## Contact

Dvir Ross
- dvirross@shenkar.ac.il
- RossDv@sce.ac.il
