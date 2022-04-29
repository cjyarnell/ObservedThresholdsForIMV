## Descriptive Analyses of PAHRC Cohort

setwd("C:/Git/PAHRC/imvObservedCriteria")
library(tidyverse)
library(ggpubr)
library(tableone)

# colours
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_blue = "cornflower blue"


###################################

# Eligibility flowchart info

# MIMIC

library(arrow)

eligibility <- read_parquet(
  file = "PAHRC_OC_Full_Eligibility"
)

nrow(eligibility) #76540
nrow(filter(eligibility,
            eligibility_admission == 0)) #9878
nrow(filter(eligibility,
            eligibility_outcome == 0)) 
nrow(filter(eligibility, 
            goc_eligible == 0)) # 106
#661 + 5 for the patients who died or were discharged before admission... 
# stay_ids 37211828, 34812932, 36359523, 36250432, 31453779

eligibility <- filter(eligibility, 
                      !(stay_id %in% c(37211828, 
                                       34812932, 
                                       36359523, 
                                       36250432, 
                                       31453779)))

nrow(filter(eligibility,
            eligibility_age == 0)) #0

second_eligibility <- filter(eligibility, 
                             eligibility_outcome == 1,
                             eligibility_admission == 1,
                             goc_eligible == 1)
nrow(second_eligibility) #65898

nrow(filter(second_eligibility,
            eligibility_trach == 0)) #1261

nrow(filter(second_eligibility, 
            ineligible_imv == 1 & eligible == 0, 
            eligibility_trach == 1)) # 22695

third_eligibility <- filter(second_eligibility,
                            ineligible_imv == 0 | eligible == 1,
                            eligibility_trach == 1)

nrow(third_eligibility) # 41942


nrow(filter(third_eligibility, fio2_geq_40 == 0)) # 33438
nrow(filter(third_eligibility, eligible_o2_device == 0)) # 36132

fourth_eligibility <- filter(third_eligibility,
                             fio2_geq_40 == 1,
                             eligible_o2_device == 1)

nrow(fourth_eligibility) # 4424

nrow(filter(fourth_eligibility, repeat_eligible_admission == 1, eligible ==0)) #230
length(unique(filter(fourth_eligibility, repeat_eligible_admission == 1, eligible ==0)$subject_id)) # 191
nrow(filter(eligibility, eligible == 1)) # 3365

eligibility2 <- read_parquet(
  file = "PAHRC_OC_Eligibility"
)


## AMDS

amd_elig <- read_parquet(
  file = "AMDS_OC_Eligibility"
)

amd_elig2 <- read_parquet(
  file = "AMDS_OC_admissions"
) %>% 
  left_join(amd_elig, by = "admissionid")

# filter the following
# 16642 17126 - IMV at moment of eligibility
# 8830 13788 - no fio2 measurements after filtering out obs concurrent with imv
# 7708 - fio2 not geq 40 after filtering 
# 18674 no measurements at eligibility or later

amd_elig <- amd_elig %>%
  filter(!(admissionid %in% c(16642, 17126, 
                              8830, 13788,
                              7708, 18674)))

23106-nrow(amd_elig) #6,049 never used NRB / NIV / IMV with FiO2 >= 40

nrow(amd_elig)

amd_elig %>% filter(operative == 1) %>% count() 

amd_elig %>% 
  filter(operative == 0) %>%
  count()

amd_elig %>% 
  filter(operative == 0) %>%
  filter(o2_device == 20) %>% 
  count() # 11105 + 4506 = 

amd_elig %>% 
  filter(operative == 0) %>%
  filter(o2_device != 20) %>% 
  count()

amd_elig %>% 
  filter(operative == 0) %>%
  filter(o2_device != 20) %>% 
  filter(goc == 3) %>%
  count()

amd_elig %>% 
  filter(operative == 0) %>%
  filter(o2_device != 20) %>% 
  filter(goc != 3 | is.na(goc)) %>%
  count()
  
amd_elig %>% 
  filter(operative == 0) %>%
  filter(o2_device != 20) %>% 
  filter(goc != 3 | is.na(goc)) %>%
  filter(eligibletime >= 1440) %>%
  count()


amd_elig %>% 
  filter(operative == 0) %>%
  filter(o2_device != 20) %>% 
  filter(goc != 3 | is.na(goc)) %>%
  filter(eligibletime < 1440) %>%
  count()


###################################

# Baseline data
pdx_t1 <- readRDS("oc_pdx_b.rds") %>%
  mutate(o2_device = ifelse(hfnc_bl == 1, "HFNC", 
                            ifelse(niv_bl == 1, "NIV",
                                   "NRB"))) %>%
  mutate(o2_device = factor(o2_device)) %>%
  rename(sex = gender) %>%
mutate(
  agegroup = cut(anchor_age, breaks = c(18, 40, 50, 60, 70, 80,100),
                 right = F, ordered_result = T,
                 labels = c("18-39",
                            "40-49",
                            "50-59",
                            "60-69",
                            "70-79",
                            "80 or more")),
  careunit = ifelse(careunit %in% c("Coronary Care Unit (CCU)",
                                         "Cardiac Vascular Intensive Care Unit (CVICU)"),
                         "Cardiac",
                         ifelse(careunit %in% c("Neuro Intermediate",
                                                "Neuro Stepdown",
                                                "Neuro Surgical Intensive Care Unit (Neuro SICU)",
                                                "Trauma SICU (TSICU)"),
                                "Neuro-trauma",
                                "Medical-surgical")),
  ethnicity = ifelse(ethnicity == "WHITE", "White",
                   ifelse(ethnicity == "BLACK/AFRICAN AMERICAN", "Black",
                          ifelse(ethnicity == "ASIAN", "Asian",
                                 ifelse(ethnicity == "HISPANIC/LATINO", "Hispanic",
                                        ifelse(ethnicity %in% c("AMERICAN INDIAN/ALASKA NATIVE",
                                                                "OTHER"), "Other", 
                                               "Unknown")))))) %>%
  mutate(death28 = (time_to_death < (24*60*28)) %in% TRUE)

table1vars <- c("agegroup",
                "sex",
                "ethnicity",
                
                "anchor_year_group",
                "careunit",
                "copd",
                "chf",
                
                "spo2_bl",
                "fio2_bl",
                "resp_rate_bl",

                "imv_ever",
                "icu_dc",
                "death_obs",
                "death28"
)

catvars <- c("agegroup", "sex","ethnicity","anchor_year_group",
             "careunit","copd","chf",
             "imv_ever","icu_dc",
             "death_obs","death28")

nonNormal <- c("spo2_bl","fio2_bl",
               "resp_rate_bl")

pdx_T1 <- CreateTableOne(vars = table1vars, data = pdx_t1,
                     factorVars = catvars,
                     includeNA = T, strata = "o2_device",
                     addOverall = T, test = F)

pdx_printT1 <- print(pdx_T1, nonnormal = nonNormal, quote = F, noSpaces = TRUE, printToggle = F)

write.csv(pdx_printT1, file = "pdx_T1.csv")

readRDS("oc_pdx_tcs_stan.rds") %>%
  filter(time >= 0 ) %>%
  filter(var == "niv") %>%
  group_by(stay_id) %>%
  summarise(niv_ever = max(value)) %>%
  summarise(sum(niv_ever),
            mean(niv_ever))

readRDS("oc_pdx_tcs_stan.rds") %>%
  filter(time >= 0 ) %>%
  filter(variable == "hfnc") %>%
  group_by(stay_id) %>%
  summarise(hfnc_ever = max(value)) %>%
  summarise(sum(hfnc_ever),
            mean(hfnc_ever))

readRDS("oc_pdx_tcs_stan.rds") %>%
  filter(time >= 0 ) %>%
  filter(variable == "stdo2") %>%
  group_by(stay_id) %>%
  summarise(stdo2_ever = max(value)) %>%
  summarise(sum(stdo2_ever),
            mean(stdo2_ever))


####

# for AMDS

amd <- readRDS("oc_amd.rds")
amd_b <- readRDS("oc_amd_b.rds")

amd %>% 
  group_by(admissionid) %>%
  filter(var == "niv") %>%
  summarise(niv_ever = max(value)) %>%
  summarise(sum(niv_ever),
            mean(niv_ever))

amd %>% 
  group_by(admissionid) %>%
  filter(var == "stdo2") %>%
  summarise(nrb_ever = max(value)) %>%
  summarise(sum(nrb_ever),
            mean(nrb_ever))

niv_bl<-
  amd %>% 
  group_by(admissionid) %>%
  filter(var == "niv", chartminute <= 0, !is.na(value)) %>%
  summarise(niv_bl = last(value))

stdo2_bl<-
  amd %>% 
  group_by(admissionid) %>%
  filter(var == "stdo2", chartminute <= 0, !is.na(value)) %>%
  summarise(stdo2_bl = last(value))
  
rr_bl<-
  amd %>% 
  group_by(admissionid) %>%
  filter(var == "resp_rate", chartminute <= 0, !is.na(value)) %>%
  summarise(resp_rate_bl = last(value))

fio2_bl<-
  amd %>% 
  group_by(admissionid) %>%
  filter(var == "fio2", chartminute <= 0, !is.na(value)) %>%
  summarise(fio2_bl = last(value))

spo2_bl<-
  amd %>% 
  group_by(admissionid) %>%
  filter(var == "spo2", chartminute <= 0, !is.na(value)) %>%
  summarise(spo2_bl = last(value))
    
amd_t1 <- 
  stdo2_bl %>%
  left_join(rr_bl, by = "admissionid") %>%
  left_join(fio2_bl, by = "admissionid") %>%
  left_join(spo2_bl, by = "admissionid") %>%
  mutate(o2_device = ifelse(stdo2_bl == 1, "NRB",
                            "NIV")) %>%
  left_join(amd_b, by = "admissionid") %>%
  select(-stdo2_bl, -id, -patientid) %>%
  rename(sex = gender,
         o2_device = o2_device.x) %>%
  mutate(sex = ifelse(sex == "Vrouw", 1, 
                      ifelse(sex == "Man", 0,NA)),
         imv = (imv_time < 7200) %in% TRUE,
         death28 = (dateofdeath < (60*24*28)) %in% TRUE,
         death_obs = (dateofdeath < 7200) & is.na(imv_time)) %>%
  mutate(icu_dc = ifelse(imv == F & 
                           ((dischargedat/60000 < 7200) %in% TRUE),
                         1, 0))

table1vars <- c("agegroup",
                "sex",
                
                "admissionyeargroup",
                "location",
                
                "spo2_bl",
                "fio2_bl",
                "resp_rate_bl",
                "goc",
                
                "imv",
                "icu_dc",
                "death_obs",
                "death28"
)

catvars <- c("agegroup", "sex","admissionyeargroup",
             "location","goc",
             "imv","icu_dc",
             "death_obs","death28")

nonNormal <- c("spo2_bl","fio2_bl",
               "resp_rate_bl")

T1 <- CreateTableOne(vars = table1vars, data = amd_t1,
                     factorVars = catvars,
                     includeNA = T, strata = "o2_device",
                     addOverall = T, test = F)

printT1 <- print(T1, nonnormal = nonNormal, quote = F, noSpaces = TRUE, printToggle = F)

write.csv(printT1, file = "amd_T1.csv")


# how many of the patients not intubated go on to be intubated later?

# Missing data at baseline


###################################

# Timevarying data

pdx <- readRDS("oc_pdx.rds")
pdx_tcs_stan <- readRDS("oc_pdx_tcs_stan.rds")


# death imv and discharge rate over time, MIMIC

lastobs <- group_by(pdx_tcs_stan, stay_id) %>%
  summarise(lastobs = max(time))

plotdf <- pdx_b %>%
  select(stay_id, death_obs_time, icu_dc_time, imv_time) %>%
  left_join(lastobs, by = "stay_id") %>%
  group_by(stay_id) %>%
  mutate(time = min(c(death_obs_time, icu_dc_time, imv_time, lastobs), na.rm = T),
         event = ifelse(!is.na(imv_time), 2, 
                        ifelse(!is.na(icu_dc_time), 3,
                               ifelse(!is.na(death_obs_time),4, 1)))) %>%
  select(-(death_obs_time:imv_time), -lastobs) %>%
  mutate(event = factor(event, levels = 1:4, 
                        labels = c("None", "IMV","DC","Death")))

library(survival)
library(survminer)
sfit_imv <- survfit(Surv(time, event=="IMV")~1, data=plotdf)
sfit_dc <- survfit(Surv(time, event=="DC")~1, data=plotdf)
sfit_death <- survfit(Surv(time, event=="Death")~1, data=plotdf)

plotdf<-  data.frame(imv = cumsum(sfit_imv$n.event),
                     dc = cumsum(sfit_dc$n.event),
                     death = cumsum(sfit_death$n.event),
                     time = sfit_imv$time) %>%
  mutate(inICU = nrow(pdx_b) - imv-dc-death) %>%
  pivot_longer(c("imv","dc","death","inICU"), names_to = "State") %>%
  mutate(State = factor(State, levels = c("death","dc","imv","inICU"),
                        labels = c("Death","ICU discharge", "IMV","ICU, no IMV"),
                        ordered = T))

plot <- ggplot(data =plotdf,
               aes(x = time/60, y = value/nrow(pdx_b), fill = State)) +
  geom_area(position = "stack") +
  theme_minimal() + 
  scale_fill_manual(values = c(c_dark, c_mid, c_light, "grey95")) +
  labs(y = "Proportion", x = "Hour since eligibility",
       title = "Patients in each state over time") +
  scale_x_continuous(limits = c(0,120),
                     breaks = seq(from = 0, to = 120, by = 24)) +
  theme(panel.grid = element_blank())

ggsave(filename = "StateTime.svg", plot = plot, height = 5, width = 5)
  