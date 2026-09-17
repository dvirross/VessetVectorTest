This Najmabadi_README20230824.txt file was generated on 20230824 by [Shahpar Najmabadi]

-------------------
GENERAL INFORMATION
-------------------


1. Title of Dataset 
Menstrual Cycles Length of Women in the USA and Canada, 1990-2013


2. Author Information

  Principal Investigator Contact Information
           Name: Joseph Stanford
           Institution: University of Utah
           Address: Department of Family & Preventive Medicine, 375 Chipeta Way, Ste. A, Salt Lake City, UT 84108
           Email: joseph.stanford@utah.edu

  Associate or Co-investigator Contact Information
           Name: Shahpar Najmabadi
           Institution: University of Utah
           Address: Department of Family & Preventive Medicine, 375 Chipeta Way, Ste. A, Salt Lake City, UT 84108
           Email: s.najmabadi@utah.edu

3. Date of data collection (single date, range, approximate date) <suggested format YYYYMMDD>
1990-2013

4. Geographic location of data collection (where was data collected?): 
Creighton Model FertilityCare Service Centers located in New Hampshire, Nebraska, Utah, New Mexico, Missouri, 
Kansas, South Carolina, Minnesota, West Virginia, California, Ohio, Texas, Pennsylvania, and Ontario, Canada.


5. Information about funding sources that supported the collection of the data:
The Robert Wood Johnson Foundation for Creighton Model Fecundity, grant number RWJF#029258. 
The Eunice Kennedy Shriver National Institute of Child Health and Human Development for Time to Pregnancy in Normal Fertility, grant number 1K23 HD0147901-01A1.
The Office of Family Planning, Office of Population Affairs, Health and Human Services, grant number 1FPRPA006035.


--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 


1. Licenses/restrictions placed on the data:
CC BY NC - Allows others to use and share your data non-commercially and with attribution.

2. Links to publications that cite or use the data:
https://pubmed.ncbi.nlm.nih.gov/32104920/
https://pubmed.ncbi.nlm.nih.gov/33990841/
https://pubmed.ncbi.nlm.nih.gov/36186844/


3. Links to other publicly accessible locations of the data:
NA

4. Links/relationships to ancillary data sets:
NA

5. Was data derived from another source? 
           If yes, list source(s):
“Creighton Model Effectiveness, Intentions, and Behaviours Assessment” (CEIBA) (2009–2013);
“Time to Pregnancy in Normal Fertility” (TTP) (2003–2006);
“Creighton Model MultiCenter Fecundability Study” (CMFS) (1990–1996).


6. Recommended citation for the data:
Menstrual Cycles Length of Women in the USA and Canada, 1990-2013. The Hive: University of Utah Research Data Repository. https://doi.org/10.7278/S50d-4gxs-s4hj



---------------------
DATA & FILE OVERVIEW
---------------------


1. File List
   A. Filename:   CrMcyclelength_share.csv      
      Short description:        
The purpose of this derived dataset was to analyze menstrual cycle lengths in relation to lunar calendar.
This datafile of start and end date of 3324 menstrual cycles of 581 women is part of a
combined dataset of three cohorts of heterosexually active women who received instruction
in the Creighton Model FertilityCare System (CrM) through centres across the United States
and Canada. The CrM has standardised protocols for teaching women how to observe, record,
and interpret daily vaginal discharge from bleeding and cervical fluid on a daily diary,
called a CrM chart, and to use these standardised observations to identify the estimated time
of ovulation and days when intercourse is likely to result in pregnancy. The cohorts included:
“Creighton Model Effectiveness, Intentions, and Behaviours Assessment” (CEIBA) (2009–2013), a 
prospective cohort of women without known subfertility, aimed to evaluate and classify pregnancy
rates and pregnancy intentions during use of the CrM; “Creighton Model MultiCenter
Fecundability Study” (CMFS) (1990–1996), a retrospective cohort of presumably fertile
and subfertile women using CrM, aimed to assess the relationship between vulvar mucus
observations and the day and cycle-specific probabilities of conception; and “Time to
Pregnancy in Normal Fertility” (TTP) (2003–2006), a parallel-randomised trial, which
aimed to assess the impact of CrM use on time to pregnancy in couples of proven
fertility trying to conceive.
Each of the cohorts aimed to include heterosexually active couples with normal fecundity.
Eligibility criteria were assessed by women’s responses to the CrM general intake form
and/or a screening questionnaire. Eligibility requirements in the original studies included
women, age 18–40 years old (upper limit of 35 years for TTP), not pregnant at entry,
having regular menstrual bleeding, and not breast feeding (CMFS and TTP), or if breast
feeding, not doing so exclusively (CEIBA). Recent users of oral contraceptives had to have
at least one menstrual bleed (CEIBA) or two menstrual bleeds (TTP) since stopping the
oral contraceptives; however, for CMFS there was no restriction for time since discontinuing
oral contraceptives. All studies also required normal menstrual patterns since last use
of depo-medroxy-progesterone acetate or a hormonal intra-uterine device. 



2. Relationship between files:   NA     




3. Additional related data collected that was not included in the current data package:  
See above publications



4. Are there multiple versions of the dataset? No 



--------------------------
METHODOLOGICAL INFORMATION
--------------------------


1. Description of methods used for collection/generation of data: 
581 women (3,324 cycles) with no known subfertility (18-40 years of age) were followed for up to 1 year. Women recorded vaginal bleeding and mucus discharge daily. We used the peak day of cervical mucus as the estimated day of ovulation and the last day of the follicular phase. We used generalised linear mixed models stratified by age and parity to describe menstrual cycle parameters.


2. Methods for processing the data: 
Data were exported from Microsoft Access, cleaned in SAS, and exported to spreadsheet format.


3. Instrument- or software-specific information needed to interpret the data:
NA

4. Standards and calibration information, if appropriate:
NA

5. Environmental/experimental conditions:
NA

6. Describe any quality-assurance procedures performed on the data:
NA

7. People involved with sample collection, processing, analysis and/or submission:
Shahpar Najmabadi
Research staff and faculty at the University of Utah, Department of Family and Preventive Medicine: Becky Crockett, Julie Fryer, Sebrena Kruletz, Iris Musso, Jared Hansen, Crystal Xu, Charles Wu, Esther Chang, Michael Lowe, Nirupma Singh, Daisy Krakowiak, Kaitlin Carruth



-----------------------------------------
DATA-SPECIFIC INFORMATION FOR: [FILENAME]
-----------------------------------------
<create sections for each dataset (or file if appropriate) included>


1. Number of variables:
Seven (7) variables

2. Number of cases/rows: 
3324.

3. Variable List
    A. Name: [new_id]
       Description: [sequential unique woman ID #]
                    Value labels if appropriate: Numeric

    B. Name: [age]
       Description: [woman's age in years (at time of first cycle)]
                    Value labels if appropriate: Numeric

    C. Name: [cycle_number]
       Description: [sequential cycle number in dataset]
                    Value labels if appropriate: Numeric

    D. Name: [cycle_start_date]
       Description: [first day of menstrual cycle]
                    Value labels if appropriate: Date

    E. Name: [cycle_end_date]
       Description: [date of last day of cycle]
                    Value labels if appropriate: Date

    F. Name: [cycle_length]
       Description: [cycle length in days (blank if conception cycle)]
                    Value labels if appropriate: Numeric, blank cells represent missing due to conception cycle

    G. Name: [conception_cycle]
       Description: [whether pregnancy occurred in this cycle, yes or no]
                    Value labels if appropriate: Yes, No, Missing

4. Missing data codes:
In The uploaded CSV file, variable cycle_length has 180 blank cells, representing missing, all related to "YES" value of the variable conception_Cycle.

5. Specialized formats of other abbreviations used
None
