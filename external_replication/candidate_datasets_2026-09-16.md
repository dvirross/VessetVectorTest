# Candidate public longitudinal cycle-length datasets for an external replication

Date: 2026-09-16. Context: the file `cycle_lengths.npz` in iurteaga/menstrual_cycle_analysis
is simulated (see `clue_audit_2026-09-16.md`); a replacement source is needed.

Requirements for the research question (rule-defined pattern counts under H_G / H_iid / H_W):
individual-level records; per-woman **consecutive** cycle lengths in integer days with
verifiable order and gap information; at least 5 cycles per woman for a useful share of women;
reproductive-age women not on hormonal contraception; provenance independent of the Fehring
(Marquette) trial; downloadable or obtainable under a documented procedure.

Network note: hive.utah.edu, stat.unipd.it, datarepository.stat.unipd.it, icpsr.umich.edu,
physionet.org, zenodo.org, dryad, osf.io, kaggle.com and all journal sites are blocked from this
sandbox. Facts below come from search-engine snippets of those pages and from the associated
papers; each must be confirmed on download.

## Tier 1 — real, individual-level, consecutive cycles, independent of Fehring

### 1. "Menstrual Cycles Length of Women in the USA and Canada, 1990–2013" (University of Utah, Hive)
- DOI 10.7278/S50d-4gxs-s4hj; creators Joseph B. Stanford and Shahpar Najmabadi; uploaded
  2023-08-24, modified 2023-11-29; licence CC BY-NC; direct download.
- Content: start and end dates of 3,324 menstrual cycles from 581 women (mean 5.7 cycles per
  woman; follow-up up to one year), pooled from three Creighton Model FertilityCare cohorts:
  CMFS 1990–1996 (retrospective), TTP 2003–2006 (trial), CEIBA 2009–2013 (prospective).
  Eligibility in the source studies: age 18–40 (≤35 in TTP), not pregnant, regular menstrual
  bleeding, not exclusively breastfeeding, no known subfertility (Najmabadi et al. 2020,
  Paediatr Perinat Epidemiol 34:318; Najmabadi et al. 2022, Hum Reprod Open hoac039).
- Used as the North-American cohort of Ecochard et al. 2024 (Sci Adv 10:eadg9646) together
  with the Marquette archive; Ecochard kept "one series of consecutive menstrual cycles" per
  woman (562 women, 3,137 cycles), so the raw file contains some non-consecutive series and
  gap handling must be pre-specified.
- Fit: closest analogue to Fehring (NFP-charting cohort, daily bleeding diary, exact start
  dates, selection towards regular cycles) with different method, investigators, sites and
  decades. Main limitation: short follow-up, so the ≥5-cycle rule will remove a large share of
  women; the exact retained n is unknown until the file is inspected.

### 2. Colombo "London" database, Data Repository of the Department of Statistical Sciences, Padua (item 19)
- DOI 10.25430/datarepository-statisticalsciences_19; depositor Bernardo Colombo; data
  described as from 1993; licence CC BY 4.0; some files require an access request through
  the repository.
- Files: `london` (all cycles of length ≤ 100 days, with daily basal body temperature; CSV
  9.28 MB, RData 0.83 MB, SAS), `anomalous` (cycles of women with any cycle > 100 days),
  `alletaev` (the london records without daily BBT plus the woman's history).
- Fit: long individual BBT-chart series (woman id, sequence of cycles, explicit handling of
  > 100-day cycles). Number of women and cycles not visible in snippets; likely to provide
  many women with long follow-up, which the Utah file lacks.

### 3. Colombo collection of cycle-biometry databases, Padua (application procedure)
- https://www.stat.unipd.it/ricerca/open-data/basi-di-dati-sulla-biometria-del-ciclo and
  https://www.stat.unipd.it/ricerca/basi-di-dati. Access: short research project plus an
  application form to the Department of Statistical Sciences.
- Includes at least: **Fertili** (European multicentre NFP study, 1992–1996, 7 centres in
  Italy, Germany, Belgium, France, Spain and Switzerland: 881 women, 7,017 cycles; inclusion
  required a history of 24–34-day cycles and mucus-monitoring experience — Colombo &
  Masarotto 2000, Demographic Research 3(5)); an Italian mucus/intercourse database (193
  women, 2,755 cycles); and the compiled European dataset used by Ecochard et al. 2024
  (2,303 women, 26,912 cycles). The 2026 WAVES paper (Sci Adv, 10.1126/sciadv.aeb1175) used
  753 women / 5,674 cycles from Fertili.
- Fit: large, long series, independent. Not instant: needs the application. Note the Fertili
  inclusion criterion selects on regularity (24–34 days), which must be reported as such.

## Tier 2 — real and obtainable, but population or structure differs

### 4. BIMORA (ICPSR 4452), Tremin sub-cohort 1998–2002
- DOI 10.3886/ICPSR04452.v2; public-use via a free ICPSR account. 156 Tremin women aged
  25–58; daily records of menstrual bleeding from 15 January to 14 July in each of five
  years (five ~6-month windows per woman, with gaps between windows), merged with the Tremin
  calendar-card history and urinary hormones.
- Fit: real consecutive cycles within each window, but a perimenopausal population and a
  windowed design; suitable as a contrasting sensitivity cohort, not as a like-for-like
  replication.

### 5. TREMIN Research Program on Women's Health (Penn State; application)
- Menstrual calendar cards since 1934, ~5,000 women, 58,296 menstrual records. The largest
  source of long consecutive series, but access is by application and not documented online.

### 6. App cohorts — not public
- Clue (BioWink): data-use agreement only (all Li/Urteaga papers). Kindara/Sympto (Symul et
  al. 2019): only aggregated tables in lasy/FAM-Public-Repo (verified: transition/emission
  probability tables and summary counts, no user records). Natural Cycles, Flo, Apple
  Women's Health Study, Oura: not available.

## Tier 3 — unsuitable

- **mcPHASES** (PhysioNet 1.0.0 / Sci Data 2026): 42 participants, 192 complete cycles in
  two 3-month intervals two years apart; ≈4.6 cycles per participant — below the ≥5 rule.
- **IEEE DataPort "Menstrual Cycle Phase Prediction"**: 42 individuals over three months.
- **Kaggle**: "Menstrual Cycle Data" (nikitabisht) is a re-upload of the Marquette archive
  (not independent); "Menstrual cycle data with factors" is declared synthetic; "Menstural
  Cycle length Dataset" has no documented provenance.
- **iurteaga/menstrual_cycle_analysis `cycle_lengths.npz`**: simulated (audit file).
- **lasy/HiddenSemiMarkov**: ships only a simulated `simple_model` example.
- **Vollman (1977)**: no digitised individual records located.

## Recommendation

Primary replication: **Utah Hive dataset (1)** — open licence, immediate download, same
construct (bleeding-start-to-bleeding-start cycle length), independent NFP cohort. Secondary,
for long follow-up: **Padua item 19 (2)**, with the **Colombo collection (3)** requested in
parallel for a larger European replication.

Pre-specification to fix *before* the Utah file is opened for pattern counting:
1. Cycle length L = start of next recorded cycle − start of this cycle (in days); if the file
   gives end dates, a cycle is "consecutive" with the next when next start = end + 1 (or the
   repository's own consecutive flag, if present).
2. Split each woman's record at any gap; keep only runs of ≥ 5 consecutive cycles; if a woman
   has several such runs, keep the longest (ties: the earliest), so each woman contributes
   one contiguous series as in Fehring.
3. No cycle-length cutoff in the primary analysis; the harmonised 18–54-day sensitivity and
   the truncation-to-first-12 sensitivity as already coded.
4. Same B = 50,000, seed 17, same engine (`scripts/patterns.py`); outputs to
   `results/external/utah/`.

Download must be done outside this sandbox (blocked hosts); place the files under
`data/external/utah/` with their SHA-256 recorded, and the pipeline will read them from there.

## Sources consulted (search-engine snippets; pages themselves blocked from the sandbox)
- https://hive.utah.edu/ (Stanford & Najmabadi dataset, DOI 10.7278/S50d-4gxs-s4hj)
- https://datarepository.stat.unipd.it/item/19
- https://www.stat.unipd.it/ricerca/open-data/basi-di-dati-sulla-biometria-del-ciclo
- https://www.science.org/doi/10.1126/sciadv.adg9646 (Ecochard et al. 2024)
- https://www.science.org/doi/10.1126/sciadv.aeb1175 (WAVES, 2026)
- https://onlinelibrary.wiley.com/doi/10.1111/ppe.12644 (Najmabadi et al. 2020)
- https://academic.oup.com/hropen/article/2022/4/hoac039/6696224 (Najmabadi et al. 2022)
- https://www.demographic-research.org/volumes/vol3/5/3-5.pdf (Colombo & Masarotto 2000)
- https://www.icpsr.umich.edu/web/DSDR/studies/4452 (BIMORA)
- https://physionet.org/content/mcphases/1.0.0/ ; https://www.nature.com/articles/s41597-026-06805-3
- https://github.com/lasy/FAM-Public-Repo (cloned and inspected)
