These scripts are used to build the cohorts.

MIMIC 
(scripts prefixed with PAHRC, Pneumonia and Acute Hypoxemic Respiratory failure Cohort)

1) O2DeliveryFiO2Table_1 to build the fio2 and oxygen delivery device table 
from existing derived and fundamental MIMIC-IV tables
2) PAHRC_OC_Eligibility to build a table of eligible ICU stay ids
3) PAHRC_OC_Timevarying to gather the timevarying data for each eligible stay
4) PAHRC_OC_Baseline to gather the baseline info for each eligible stay

AmsterdamUMCdb
(scripts prefixed with AMDS)
1) Preprocessing scrips: AMDS_FiO2, AMDS_goc, AMDS_ph, AMDS_po2pco2, AMDS_vitals
2) AMDS_OC_Eligibility
3) AMDS_OC_timevarying

Drives and table names may need to be adjusted to ensure smooth flow.

After these scripts are run, download the tables corresponding to 2/3/4 
and 2/3 above then go to the R code in the CreateCohort folder