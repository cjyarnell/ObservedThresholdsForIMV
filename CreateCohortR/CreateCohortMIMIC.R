###################################################
#
# finalize cohort for analysis

##################################################
# MIMIC-IV aka PAHRC first
# Pneumonia Acute Hypoxemic Respiratory Failure Cohort
# for the investigation of observed criteria for invasive ventilation
# builds cohort database from 2 parquet files:
# PAHRC_OC_Timevarying_
# PAHRC_OC_Baseline

library(tidyverse)
library(caret)
library(arrow)

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

setwd("C:/Git/PAHRC/imvObservedCriteria")
pdx <- read_parquet(
  file = "PAHRC_OC_Timevarying"
)

pdx_b <- read_parquet(
  file = "PAHRC_OC_Baseline"
) 

# remove five patients who have death or discharge charted before eligibility
pdx_b <- filter(pdx_b, !(stay_id %in% c(37211828, 34812932, 36359523, 36250432, 31453779)))
pdx <- filter(pdx, !(stay_id %in% c(37211828, 34812932, 36359523, 36250432, 31453779)))

# fix 2 patients with NA for chf and NA for copd
sum(is.na(pdx_b$chf))
sum(is.na(pdx_b$copd))

pdx_b$copd[is.na(pdx_b$copd)] <- 0
pdx_b$chf[is.na(pdx_b$chf)] <- 0

# have to tidy up the death_obs_time and icu_dc_time variables
# and the death_obs and icu_dc variables
# (make sure no deaths after discharge are included and vice versa)

pdx_b <- pdx_b %>%
  mutate(death_obs = ifelse(death_obs == 1 & icu_dc == 1
                            & death_obs_time > icu_dc_time,
                            0, death_obs),
          icu_dc = ifelse(death_obs == 1 & icu_dc == 1
                         & death_obs_time <= icu_dc_time,
                         0, icu_dc)) %>%
  mutate(death_obs_time = ifelse(icu_dc == 1, NA, death_obs_time),
         icu_dc_time = ifelse(death_obs == 1, NA, icu_dc_time))


#################################################

# BASELINE

# create transformed version of the baseline data

# push admission times to hourly

pdx_b$intime <- as.numeric(format(strptime(pdx_b$intime,"%Y-%m-%d %H:%M:%S"),'%H'))


#################################################

# TIMEVARYING

# combine GCS into one variable

# convert GCS EVM to single score

gcs <- filter(pdx,
              grepl("GCS", variable)) %>%
  group_by(stay_id, time) %>%
  summarise(variable = "GCS",
            value = sum(value),
            count = n()) %>%
  # only keep the ones that have measurements for all three GCS components
  filter(count == 3) %>%
  select(-count)

# add back to main model

pdx <- pdx %>%
  bind_rows(gcs) %>%
  filter(!(variable %in% c("GCS_eyes", "GCS_motor", "GCS_verbal")))

pdx <- pdx %>%
  filter(time <= 120*60, time > -12*60)

# transformations for continuous variables to facilitate imputation

log_group <- c("albumin", 
               "alt",
               "bicarb",
               "bili",
               "cr",
               "dbp",
               "glucose",
               "hb",
               "heart_rate",
               "lactate",
               "mbp",
               "po2",
               "pco2",
               "ph",
               "plt",
               "resp_rate",
               "sbp",
               "temperature",
               "trop",
               "uo_6h_rate",
               "wbc")

transform <- function(row, log_group){
  variable <- row[3]
  value <- as.numeric(row[4])
  if(variable %in% log_group){log(value)}
  # use logit transformation functions for GCS, spo2, and fio2
  # because they have a restricted range 
  else if (variable == "GCS"){logit((value-2.9)/12.2)}
  else if (variable == "spo2"){logit(value/101)} # take the floor when transforming back
  else if (variable == "fio2"){logit((value-20.9)/79.3)} # take the floor when transforming back
  else value
}

pdx_transformed <- pdx  
pdx_transformed$value <- apply(pdx, 1, transform, log_group)

# 26 NA from log transformations, remove
pdx_transformed <- filter(pdx_transformed,
                          !is.na(value))

cont_vars <- c("albumin", 
               "alt",
               "bicarb",
               "bili",
               "cr",
               "dbp",
               "glucose",
               "hb",
               "heart_rate",
               "lactate",
               "mbp",
               "po2",
               "pco2",
               "ph",
               "plt",
               "resp_rate",
               "sbp",
               "temperature",
               "trop",
               "uo_6h_rate",
               "wbc", 
               "spo2",
               "fio2",
               "GCS")

pdx_tcs <- pdx_transformed


# the time series ends when imv has been initiated,
# remove the imv observations

pdx_tcs <- filter(pdx_tcs, !(variable == "o2_device" & value == 6))

# make sure no observations after death or discharge
pdx_tcs <- pdx_tcs %>%
  left_join(select(pdx_b, stay_id, death_obs, death_obs_time, icu_dc, icu_dc_time), by = "stay_id") %>%
  filter(
    # no discharge and no death
    (death_obs == 0 & icu_dc == 0) |
      # death   
      (death_obs == 1 & time <= death_obs_time) |
      # discharge 
      (icu_dc == 1 & time <= icu_dc_time)) %>% 
  filter(
    # remove all observations at same time as death
    !(death_obs == 1 & time == death_obs_time & variable != "death")
  ) %>%
  select(-death_obs_time, -icu_dc_time,
         -death_obs, -icu_dc) 

nrow(filter(pdx_tcs, variable == "death", value == 1)) # 230 matches outcome table
nrow(filter(pdx_tcs, variable == "icu_dc", value == 1)) # 1964 matches outcome table

lastobs <- pdx_tcs %>% 
  group_by(stay_id) %>%
  summarise(lastobs = max(time)) 

# binary o2 device version (categorical too slow...)
# the argument is that it indicates all recently used oxygen devices
# it's not perfect but since the latent variables are correlated 
# it is somewhat addressed by the data

niv <- filter(pdx_tcs, variable == "o2_device") %>%
  mutate(variable = "niv",
         value = ifelse(value == 5, 1, 0))

hfnc <- filter(pdx_tcs, variable == "o2_device") %>%
  mutate(variable = "hfnc",
         value = ifelse(value == 4, 1, 0))

stdo2 <- filter(pdx_tcs, variable == "o2_device") %>%
  mutate(variable = "stdo2",
         value = ifelse(value == 3 | 
                          value == 2, 1, 0))


# nasal prongs and room air is the reference category (0,0,0)

pdx_tcs <- pdx_tcs %>%
  filter(!(variable == "o2_device")) %>%
  bind_rows(niv) %>%
  bind_rows(hfnc) %>%
  bind_rows(stdo2)

# record the means and standard deviations of 
# transformed variables for reconstruction later
mean_sd <- filter(pdx_tcs, variable %in% cont_vars) %>%
  group_by(variable) %>%
  summarise(mean = mean(value),
            sd = sd(value))

saveRDS(mean_sd, file = "oc_mean_sd2.rds")

for (j in 1:nrow(mean_sd)){
  var <- as.character(mean_sd[j,1])
  mean <- as.double(mean_sd[j,2])
  sd <- as.double(mean_sd[j,3])
  
  pdx_tcs <- mutate(
    pdx_tcs,
    value = ifelse(variable == var,
                   (value-mean)/sd,
                   value))
}


# filter out the variables we are not using in the Stan model
# computational constraints limit us to a subset of variables
# maybe one day with variational methods we can include more
# but with full MCMC it is not feasible

pdx_tcs_stan <- filter(pdx_tcs,
                       variable %in% c("resp_rate",
                                       "spo2",
                                       "heart_rate",
                                       "sbp",
                                       "fio2",
                                       "po2",
                                       "ph",
                                       "pco2",
                                       "GCS",
                                       "wob",
                                       "pressor",
                                       "hfnc",
                                       "niv",
                                       "stdo2"),
                       time > -12*60)



pdx_tcs_stan <- 
  pivot_wider(pdx_tcs_stan,
              names_from = variable,
              values_from = value,
              values_fill = NA) %>%
  group_by(stay_id) %>%
  fill(fio2, pressor, hfnc, niv, stdo2) %>%
  pivot_longer(cols = sbp:stdo2,
               names_to = "var",
               values_to = "value")

#########################################

# Split cohort into manageable bites
# to find Gaussian process hyperparameters
# take a random contiguous 48h from every patient

K_total <- nrow(pdx_b)
set.seed(20210205)
interval_length <- 48*60
K_subset <- sample(1:K_total, size = 400, replace = F)
patients <- data.frame(id = 1:K_total, 
                       stay_id = pdx_b$stay_id)

patients2 <- filter(patients, id %in% K_subset) %>%
  mutate(id = 1:length(K_subset))

pdx_b <- left_join(pdx_b, patients, by = "stay_id")

library(cmdstanr)

# prep data for stan model

# data for interpolating each patient

cvars <- c("sbp", "resp_rate",  "heart_rate", "spo2","fio2", 
           "GCS", "po2", "ph","pco2") # continuous variables 
bvars <- c("wob","pressor", "hfnc","niv","stdo2")

D_c <- length(cvars)
D_b <- length(bvars)
D = D_c + D_b

dft <- left_join(patients, pdx_tcs_stan, by = "stay_id") %>%
  mutate(var = factor(var, levels = c(cvars,bvars), labels = 1:D)) %>%
  mutate(time = time/(24*60)) %>%
  arrange(id, var, time) %>%
  mutate(value = ifelse(is.na(value) & var == 11, 0, value),
         var = as.integer(as.character(var))) 

dft_pred <- filter(dft, is.na(value)) %>%
  filter(var %in% c(1:4,6,10)) %>%
  filter(time >= 0)

dft <- filter(dft, !is.na(value))

k_obs <- dft$id

# number of observations per patient
N_k <- dft %>%
  group_by(id) %>%
  summarise(N_k = n()) %>%
  select(id,N_k)

N_kc <- dft %>%
  filter(var <= D_c) %>%
  group_by(id) %>%
  summarise(N_kc = n()) %>%
  select(id,N_kc)

N_kb <- dft %>%
  filter(var > D_c) %>%
  group_by(id) %>%
  summarise(N_kb = n()) %>%
  select(id,N_kb)

Ns <- left_join(N_k, N_kc, by = "id") %>%
  left_join(N_kb, by = "id")

d_obs <- dft$var
d_obs_matrix <- matrix(0,nrow = nrow(dft), ncol = D)
for (n in 1:nrow(dft)){d_obs_matrix[n,d_obs[n]] = 1}

y_c <- dft$value[d_obs <= D_c]
y_b <- dft$value[d_obs > D_c]

x_pred <- dft_pred$time 

N_pred <- nrow(dft_pred)

d_pred <- dft_pred$var
#d_pred_matrix <- matrix(0,nrow = nrow(dft_pred), ncol = D)
#for (n in 1:nrow(dft_pred)){d_pred_matrix[n,d_pred[n]] = 1}

N_k_pred <- dft_pred %>%
  group_by(id) %>%
  summarise(n = n()) %>%
  select(n)
N_k_pred <- as.integer(unlist(N_k_pred))

interp_dat <- list(N_obs = nrow(dft),
                   K = length(unique(dft$id)),
                   D = D,
                   N_k = Ns$N_k,
                   N_kc = Ns$N_kc,
                   N_kb = Ns$N_kb,
                   
                   x_obs = dft$time,
                   d_obs = dft$var,
                   d_obs_matrix = d_obs_matrix,
                   k_obs = dft$id,
                   
                   D_c = D_c,
                   N_c = nrow(filter(dft, var <= D_c)),
                   y_c = y_c,
                   
                   D_b = D_b,
                   N_b = nrow(filter(dft, var > D_c)),
                   y_b = y_b,
                   
                   N_pred = N_pred,
                   k_pred = dft_pred$id,
                   N_k_pred = N_k_pred,
                   x_pred = x_pred,
                   d_pred = d_pred,

                   patients = patients
)

saveRDS(interp_dat, file = "C:/Users/chris/OneDrive/Documents/LargeDataFiles/interpolation_dat_obsonly.rds")

saveRDS(pdx_b, "oc_pdx_b.rds")
saveRDS(pdx, "oc_pdx.rds")
saveRDS(pdx_tcs, "oc_pdx_tcs.rds")
saveRDS(pdx_tcs_stan, "oc_pdx_tcs_stan.rds")

# version to fit initial GP

dft <- left_join(patients2, select(dft,-id), by = "stay_id") %>%
  filter(time < interval_length,
         time >= 0)

# prep data for stan model

M = 15
L = 1.5

k_obs <- dft$id

# number of observations per patient
N_k <- dft %>%
  group_by(id) %>%
  summarise(N_k = n()) %>%
  select(id,N_k)

N_kc <- dft %>%
  filter(var <= D_c) %>%
  group_by(id) %>%
  summarise(N_kc = n()) %>%
  select(id,N_kc)

N_kb <- dft %>%
  filter(var > D_c) %>%
  group_by(id) %>%
  summarise(N_kb = n()) %>%
  select(id,N_kb)

Ns <- left_join(N_k, N_kc, by = "id") %>%
  left_join(N_kb, by = "id")

d_obs <- dft$var
d_obs_matrix <- matrix(0,nrow = nrow(dft), ncol = D)
for (n in 1:nrow(dft)){d_obs_matrix[n,d_obs[n]] = 1}

y_c <- dft$value[d_obs <= D_c]
y_b <- dft$value[d_obs > D_c]

interp_GP_dat <- list(N_obs = nrow(dft),
                      K = length(unique(dft$id)),
                      D = D,
                      N_k = Ns$N_k,
                      N_kc = Ns$N_kc,
                      N_kb = Ns$N_kb,
                      
                      x_obs = dft$time,
                      d_obs = dft$var,
                      d_obs_matrix = d_obs_matrix,
                      k_obs = dft$id,
                      
                      D_c = D_c,
                      N_c = nrow(filter(dft, var <= D_c)),
                      y_c = y_c,
                      
                      D_b = D_b,
                      N_b = nrow(filter(dft, var > D_c)),
                      y_b = y_b,
                      
                      M = M,
                      L = L
)

write_stan_json(data = interp_GP_dat, file = "C:/Users/chris/OneDrive/Documents/LargeDataFiles/data_for_GP_interpolation_obsonly.json")


# missing data table

pdx_tcs_stan %>%   group_by(var) %>%
  summarise(nasum = sum(is.na(value)),
            namean = mean(is.na(value)),
            count = sum(!is.na(value))) %>%
  write_csv(file = "pdx_missing.csv")

lastobs <- pdx_tcs %>% 
  group_by(stay_id) %>%
  summarise(lastobs = max(time)) 

pdx_tcs_stan %>% filter(!is.na(value)) %>%
  group_by(stay_id, var) %>%
  summarise(count = n()) %>%
  left_join(lastobs, by = "stay_id") %>%
  mutate(lastobs = lastobs+1) %>%
  group_by(var) %>%
  summarise(mean(lastobs)/mean(count)) %>%
  write_csv(file = "pdx_time_btw_obs.csv")

# fixed forward fill imputation dataset

pdx_fi <- 
  pdx_tcs_stan %>%  
  pivot_wider(names_from = "var",
              values_from = "value") %>%
  fill(c(sbp:pressor, GCS:stdo2)) %>%
  mutate(pressor = ifelse(is.na(pressor),0,pressor)) %>%
  pivot_longer(cols = sbp:stdo2,
               names_to = "variable") %>%
  filter(time >= 0, !is.na(value)) %>%
  left_join(mean_sd, by = "variable") %>%
  rename(var = variable) %>%
  mutate(value = ifelse(var %in% cvars, value*sd + mean, value)) %>%
  mutate(value = ifelse(var %in% log_group,
                        exp(value),
                        value)) %>%
  mutate(value = ifelse(var == "spo2",
                        round((inv_logit(value))*101),
                        value)) %>%
  mutate(value = ifelse(var == "GCS",
                        round((inv_logit(value))*12.2+2.9),
                        value)) %>%
  mutate(value = ifelse(var == "fio2",
                        round((inv_logit(value))*79.3+20.9),
                        value)) %>%
  select(-mean, -sd) %>%
  mutate(time = time/60)

saveRDS(pdx_fi, file = "oc_pdx_fi.rds")

