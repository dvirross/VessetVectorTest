# External dataset: University of Utah cycle-length deposit

**Source**: Stanford, J. B., & Najmabadi, S. (2023). *Menstrual Cycles Length of Women in the
USA and Canada, 1990–2013* [Data set]. The Hive, University of Utah.
https://doi.org/10.7278/S50d-4gxs-s4hj — licence CC BY-NC.

Three Creighton Model FertilityCare cohorts (CMFS 1990–1996, TTP 2003–2006, CEIBA 2009–2013):
start and end dates of 3,324 menstrual cycles from 581 women aged 18–40 with regular bleeding and
no known subfertility (Najmabadi et al. 2020, *Paediatr Perinat Epidemiol* 34:318–327;
Najmabadi et al. 2022, *Hum Reprod Open* 2022(4):hoac039; Ecochard et al. 2024, *Sci Adv*
10:eadg9646).

hive.utah.edu is not reachable from the analysis sandbox. To run the replication:

1. Download the data file(s) from the DOI above and place the original file(s) in this
   directory (`external_replication/data/raw/`).
2. `python external_replication/scripts/prepare_utah.py external_replication/data/raw/<file>` → `sequences.csv`
   and `preprocessing.json` (rules are pre-specified in the script docstring).
3. `python external_replication/scripts/run_external.py utah` → `external_replication/results/` (replicates, summary,
   validation, figures, comparison with Fehring).

`sequences.csv` has the analysed columns of `data/FilteredData.csv` (ClientID, CycleNumber,
LengthofCycle); `CycleNumber` is contiguous within every woman.
