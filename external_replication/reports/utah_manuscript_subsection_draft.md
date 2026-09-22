# Draft subsection (not yet inserted in the manuscript)

> **Superseded numbers (2026-09-22).** After this report was written the Week-Dilug rule was generalised from H = 30 to all one-weekday-progression anchors H ≡ 2 (mod 7) = {23, 30, 37, …} and the whole pipeline was rerun. Utah Week-Dilug count 40 → 42 (2 events at H = 37); joint D_M 9.99/9.53/1.20 → 10.04/9.56/1.22 (H_W p .921 → .916). Current numbers: `external_replication/results/comparison.md` and `summary.json`; the manuscript and SI Section S4 are up to date.


## External Replication in an Independent Creighton Model Cohort

To test whether the contrast between the pooled and within-woman nulls is a property of the Marquette
sample, we repeated the complete analysis, without modification, on a publicly deposited cohort
collected by different investigators with a different fertility-awareness method and in different
decades: "Menstrual Cycles Length of Women in the USA and Canada, 1990--2013" (Stanford \& Najmabadi,
2023; The Hive, University of Utah, doi:10.7278/S50d-4gxs-s4hj; CC BY-NC). The deposit gives the start
and end dates of 3,324 cycles from 581 women aged 18--40 with regular bleeding and no known
subfertility, pooled from three Creighton Model FertilityCare cohorts (CMFS 1990--1996, TTP 2003--2006,
CEIBA 2009--2013; Najmabadi et al., 2020). Preprocessing rules were fixed before any pattern was
counted: cycle length is the recorded start-to-start interval (equal to end $-$ start $+$ 1 for every
recorded cycle); the 180 conception cycles, all terminal, have no length and were excluded; a woman's
record was split wherever the next start did not follow the previous end (186 such gaps, each an
excluded cycle in the source), and gaps were never bridged; each woman contributed her longest run of
at least five consecutive cycles, as in the Marquette analysis. No cycle-length criterion was applied.
The analysed sample is 270 women and 2,139 cycles (5--15 cycles per woman, $M = 7.9$; $L$ 17--98 days,
$M = 30.0$, SD 6.2, median 29). All code, seeds ($B = 50{,}000$, seed 17) and simulation sizes were
those of the main analysis; one implementation detail was corrected beforehand (the Week anchor set is
now defined by the calendar rule $L \equiv 0 \pmod 7$ rather than derived from the observed minimum,
which leaves every Marquette result unchanged).

The observed counts are $\mathbf{O} = (41, 76, 35, 40, 41)$. Under $\HG$ all five counts exceed their
null means: Haflaga $41$ versus $12.4$ ($z = 8.33$, two-sided $p < .001$, no replicate reached 41),
Dilug $z = 3.01$, Week-Dilug $z = 3.02$, Dilug-in-Dilug $z = 2.76$ (Holm-adjusted $p$ .021--.024) and
Week $z = 1.89$ ($p = .083$); the joint statistic is $\DM = 9.99$ ($p < .001$; split-batch 10.03), and
$\HI$ gives $\DM = 9.53$. Under $\HW$ nothing remains: every $|z| \le 1.13$, all two-sided $p \ge .31$,
$\DM = 1.20$ ($p = .92$; split-batch 1.20, $p = .92$), and the exploratory three-pattern test gives
$\DM = 1.17$ ($p = .72$). The Dilug-in-Dilug deficit that was borderline in the Marquette data
($z = -2.07$, Holm $p = .19$) is absent ($z = -0.08$). The conclusion is insensitive to leaving any
woman out (joint $p$ .86--.97), to truncation at twelve cycles (2,109 cycles, $\DM = 1.34$, $p = .88$),
and to restricting the support to the 18--54-day range of the Marquette data (28 cycles removed,
$\DM = 1.15$, $p = .93$). Empirical size of all procedures is at or below nominal (marginal
.028--.048, Holm .034--.039, joint .046--.052) and the power of the $\HW$ tests against the AR(1) and
persistence alternatives is close to that in the Marquette data (joint power .44, .99 and 1.00 at
$\rho = .25, .50, .75$; .97 and 1.00 at $q = .10$ and $.20$). Of the exploratory between-woman
associations, Haflaga $\times$ Dilug-in-Dilug did not replicate ($r_s = -.09$); the only non-null pair
was Haflaga $\times$ Week ($r_s = .27$, Holm $p = .001$), which is partly a containment effect.

The independent cohort therefore reproduces, more sharply, the pattern that motivates this case study:
a count vector far outside the range of pooled exchangeability (here with all five rules elevated) that
is entirely unremarkable once each woman's own multiset of cycle lengths is conditioned on. As in the
main analysis, this is consistent with between-woman heterogeneity in cycle regularity as the source
of the population-level signal, and it says nothing about temporal ordering beyond the exact count
statistics tested.
