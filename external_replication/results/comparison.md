# utah vs Fehring — side-by-side

## Dataset characteristics

| | Fehring | utah |
|---|---|---|
| women | 118 | 270 |
| cycles | 1554 | 2139 |
| cycles/woman mean | 13.17 | 7.92 |
| cycles/woman SD | 6.06 | 2.68 |
| cycles/woman min | 5 | 5 |
| cycles/woman max | 45 | 15 |
| L mean | 29.27 | 29.98 |
| L SD | 3.86 | 6.16 |
| L median | 28.5 | 29.0 |
| L min | 18 | 17 |
| L max | 54 | 98 |
| share L >= 35 | 0.096 | 0.12 |
| Week anchors (H) | [22, 29, 36, 43, 50] | [22, 29, 36, 43, 50, 57, 64, 71, 78, 85, 92, 99] |

## Observed five-pattern counts

| Pattern | Fehring | utah |
|---|---|---|
| Haflaga | 26 | 41 |
| Dilug | 61 | 76 |
| Week | 32 | 35 |
| Week-Dilug | 18 | 42 |
| Dilug-in-Dilug | 27 | 41 |

## Marginal tests under H_G global permutation

| Pattern | Fehring μ (σ) | z | p | Holm | utah μ (σ) | z | p | Holm |
|---|---|---|---|---|---|---|---|---|
| Haflaga | 11.31 (3.27) | +4.50 | 0.0002 | 0.001 | 12.42 (3.43) | +8.33 | 0.0000 | 0.000 |
| Dilug | 48.92 (6.57) | +1.84 | 0.0850 | 0.340 | 54.95 (6.99) | +3.01 | 0.0056 | 0.017 |
| Week | 27.02 (4.01) | +1.24 | 0.2647 | 0.794 | 27.09 (4.18) | +1.89 | 0.0834 | 0.083 |
| Week-Dilug | 16.61 (3.3) | +0.42 | 0.7744 | 1.000 | 28.25 (4.26) | +3.23 | 0.0024 | 0.010 |
| Dilug-in-Dilug | 27.17 (5.0) | -0.03 | 1.0000 | 1.000 | 27.08 (5.04) | +2.76 | 0.0120 | 0.024 |

## Marginal tests under H_iid multinomial

| Pattern | Fehring μ (σ) | z | p | Holm | utah μ (σ) | z | p | Holm |
|---|---|---|---|---|---|---|---|---|
| Haflaga | 11.48 (3.31) | +4.39 | 0.0003 | 0.001 | 12.56 (3.47) | +8.20 | 0.0000 | 0.000 |
| Dilug | 48.86 (6.59) | +1.84 | 0.0820 | 0.328 | 54.84 (7.05) | +3.00 | 0.0052 | 0.021 |
| Week | 27.09 (5.0) | +0.98 | 0.3732 | 1.000 | 27.17 (5.05) | +1.55 | 0.1559 | 0.156 |
| Week-Dilug | 16.73 (3.99) | +0.32 | 0.8173 | 1.000 | 28.4 (5.16) | +2.64 | 0.0163 | 0.042 |
| Dilug-in-Dilug | 27.3 (5.07) | -0.06 | 1.0000 | 1.000 | 27.09 (5.09) | +2.73 | 0.0141 | 0.042 |

## Marginal tests under H_W within-woman permutation

| Pattern | Fehring μ (σ) | z | p | Holm | utah μ (σ) | z | p | Holm |
|---|---|---|---|---|---|---|---|---|
| Haflaga | 31.84 (4.68) | -1.25 | 0.2536 | 1.000 | 35.6 (4.78) | +1.13 | 0.3128 | 1.000 |
| Dilug | 67.04 (7.38) | -0.82 | 0.4596 | 1.000 | 75.71 (7.75) | +0.04 | 1.0000 | 1.000 |
| Week | 33.79 (3.88) | -0.46 | 0.7489 | 1.000 | 35.61 (3.97) | -0.15 | 0.9857 | 1.000 |
| Week-Dilug | 21.46 (3.3) | -1.05 | 0.3702 | 1.000 | 40.27 (4.25) | +0.41 | 0.7642 | 1.000 |
| Dilug-in-Dilug | 39.12 (5.84) | -2.07 | 0.0383 | 0.192 | 41.48 (5.95) | -0.08 | 1.0000 | 1.000 |

## Joint tests

| Null | Statistic | Fehring | utah |
|---|---|---|---|
| H_G | 5-pattern D_M | 5.074 | 10.04 |
| H_G | 5-pattern p | 0.0002 | 0.0 |
| H_G | split-batch D_M | 5.066 | 10.074 |
| H_G | split-batch p | 0.0002 | 0.0 |
| H_G | 3-pattern D_M (exploratory) | 4.499 | 9.147 |
| H_G | 3-pattern p | 0.0003 | 0.0 |
| H_G | covariance condition number | 4.57 | 4.52 |
| H_iid | 5-pattern D_M | 4.865 | 9.559 |
| H_iid | 5-pattern p | 0.0003 | 0.0 |
| H_iid | split-batch D_M | 4.868 | 9.551 |
| H_iid | split-batch p | 0.0003 | 0.0 |
| H_iid | 3-pattern D_M (exploratory) | 4.39 | 8.839 |
| H_iid | 3-pattern p | 0.0005 | 0.0 |
| H_iid | covariance condition number | 4.35 | 4.45 |
| H_W | 5-pattern D_M | 2.972 | 1.216 |
| H_W | 5-pattern p | 0.1146 | 0.9161 |
| H_W | split-batch D_M | 2.994 | 1.213 |
| H_W | split-batch p | 0.1074 | 0.9163 |
| H_W | 3-pattern D_M (exploratory) | 2.719 | 1.188 |
| H_W | 3-pattern p | 0.0591 | 0.7048 |
| H_W | covariance condition number | 5.17 | 3.96 |

## Leave-one-woman-out (H_W)

| | Fehring | utah |
|---|---|---|
| joint p range | 0.0506–0.2069 | 0.8557–0.9702 |
| # LOO joint p < .05 | 0 | 0 |
| marginal p min | [0.1206, 0.3427, 0.5338, 0.2427, 0.0142] | [0.248, 0.8409, 0.7918, 0.6117, 0.7952] |
| marginal p max | [0.3867, 0.6004, 0.9298, 0.5771, 0.073] | [0.4665, 1.0, 1.0, 0.9388, 1.0] |

## Harmonised-support sensitivity (18 <= L <= 54), external dataset

blocks 258 from 258 women, 2050 cycles; 28 cycles removed; 0 women split at a gap.

| Null | 5-pattern D_M | p | split p | marginal two-sided p (Haf, Dil, Wk, W-D, DiD) | Holm |
|---|---|---|---|---|---|
| H_G | 9.324 | 0.0 | 0.0 | [0.0, 0.0136, 0.1153, 0.0047, 0.0243] | [0.0002, 0.0407, 0.1153, 0.0189, 0.0486] |
| H_iid | 8.884 | 0.0 | 0.0 | [0.0, 0.0118, 0.2049, 0.0224, 0.028] | [0.0002, 0.047, 0.2049, 0.0671, 0.0671] |
| H_W | 1.198 | 0.9214 | 0.9232 | [0.3158, 1.0, 1.0, 0.7334, 1.0] | [1.0, 1.0, 1.0, 1.0, 1.0] |

## Simulation validation

| | Fehring | utah |
|---|---|---|
| size H_G: marginal | [0.029, 0.047, 0.034, 0.034, 0.046] | [0.039, 0.042, 0.03, 0.036, 0.046] |
| size H_G: Holm / joint plug-in / split | 0.034 / 0.050 / 0.053 | 0.041 / 0.051 / 0.052 |
| size H_iid: marginal | [0.048, 0.038, 0.04, 0.028, 0.046] | [0.028, 0.039, 0.038, 0.04, 0.048] |
| size H_iid: Holm / joint plug-in / split | 0.030 / 0.048 / 0.046 | 0.038 / 0.053 / 0.047 |
| size H_W: marginal | [0.043, 0.041, 0.043, 0.032, 0.04] | [0.045, 0.042, 0.044, 0.029, 0.035] |
| size H_W: Holm / joint plug-in / split | 0.036 / 0.052 / 0.051 | 0.032 / 0.050 / 0.051 |
| power ar1 0.0: Holm / joint | 0.025 / 0.050 | 0.028 / 0.048 |
| power ar1 0.25: Holm / joint | 0.265 / 0.525 | 0.220 / 0.455 |
| power ar1 0.5: Holm / joint | 0.890 / 0.998 | 0.885 / 0.995 |
| power ar1 0.75: Holm / joint | 1.000 / 1.000 | 1.000 / 1.000 |
| power persist 0.0: Holm / joint | 0.050 / 0.060 | 0.028 / 0.052 |
| power persist 0.1: Holm / joint | 0.938 / 0.983 | 0.945 / 0.985 |
| power persist 0.2: Holm / joint | 1.000 / 1.000 | 1.000 / 1.000 |
| power persist 0.3: Holm / joint | 1.000 / 1.000 | 1.000 / 1.000 |
| first-12 truncation: cycles, D_M, p | 1289, 2.52, 0.272 | 2109, 1.33, 0.884 |

## Exploratory pairwise (per-woman rates, Spearman, permutation p, Holm over 10 pairs)

| Pair | Fehring r / p / Holm | utah r / p / Holm |
|---|---|---|
| Haflaga x Week | 0.1 / 0.277 / 1.0 | 0.27 / 0.0001 / 0.001 |
| Haflaga x Dilug-in-Dilug | 0.27 / 0.0038 / 0.038 | -0.09 / 0.1445 / 1.0 |
| Week-Dilug x Dilug-in-Dilug | 0.07 / 0.4802 / 1.0 | 0.08 / 0.2127 / 1.0 |
| Haflaga x Week-Dilug | 0.23 / 0.0113 / 0.102 | 0.07 / 0.213 / 1.0 |
| Dilug x Dilug-in-Dilug | -0.02 / 0.819 / 1.0 | 0.07 / 0.2228 / 1.0 |
| Haflaga x Dilug | 0.04 / 0.668 / 1.0 | -0.07 / 0.2434 / 1.0 |
| Dilug x Week | 0.09 / 0.3374 / 1.0 | -0.05 / 0.368 / 1.0 |
| Week x Week-Dilug | 0.06 / 0.4755 / 1.0 | -0.02 / 0.8113 / 1.0 |
| Dilug x Week-Dilug | -0.12 / 0.2059 / 1.0 | -0.0 / 0.9497 / 1.0 |
| Week x Dilug-in-Dilug | 0.09 / 0.3262 / 1.0 | -0.0 / 0.9802 / 1.0 |
