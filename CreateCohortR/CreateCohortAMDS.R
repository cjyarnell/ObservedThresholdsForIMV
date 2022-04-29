###################################################
#
# AmsterdamUMCdb cohort construction

library(tidyverse)
library(caret)
library(arrow)

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

setwd("C:/Git/PAHRC/imvObservedCriteria")
amd <- read_parquet(
  file = "AMDS_OC_Timevarying"
)
amd_b <- read_parquet(
  file = "AMDS_OC_admissions"
)
amd_elig <- read_parquet(
  file = "AMDS_OC_Eligibility"
)

head(amd)

length(unique(amd$admissionid)) #1410

# filter charttimes > 7200 (120*60)

amd <- amd %>%
  filter(chartminute < 7200)

# filter everything after invasive ventilation

imv_ever <- amd %>%
  filter(var == "o2_device", value == 20) %>%
  group_by(admissionid) %>%
  summarise(imv_time = first(chartminute, order_by = chartminute))

length(unique(amd$admissionid))
temp <- unique(amd$admissionid)

amd <- amd %>%
  left_join(imv_ever, by = "admissionid") %>%
  filter(((chartminute < imv_time) %in% TRUE) | (is.na(imv_time))) %>%
  select(-imv_time)

length(unique(amd$admissionid))
# 2 patients were invasively ventilated at supposedly eligibility moment
# they are no longer eligible (16642 17126)

temp[!(temp %in% unique(amd$admissionid))]

# remove everything more than 24h before eligibility
amd <- filter(amd, chartminute > -1440)

# remove / fix impossible values

# fio2 between 21 and 100 keep
# between 0.21 and 1 multiply by 100
# 1 to 21 drop
amd <- amd %>%
  filter(!(var == "fio2" & value > 1 & value < 21)) %>%
  mutate(value = ifelse(var == "fio2" &
                          value <= 1 & value >= 0.21,
                        value*100,value)) %>%
  # drop mbp < 10 or > 200
  filter(!(var == "mbp" & 
             (value < 10 | value > 200))) %>%
  # drop ph < 6.5 or > 8
  filter(!(var == "ph" & 
             (value < 6.5 | value > 8))) %>%
  # fix pco2 of 0
  filter(!(var == "pco2" & value < 1)) %>%
  # remove RR > 70 or < 1
  filter(!(var == "resp_rate" & 
             (value < 1 | value > 70))) %>%
  # remove HR 0
  filter(!(var == "heart_rate" & value < 1)) %>%
  # remove spo2 > 100 or < 1
  filter(!(var == "spo2" & 
             (value > 100 | value < 1))) #%>%
#  group_by(var) %>%
#  summarise(min = min(value, na.rm = T),
#            max = max(value, na.rm = T))


# take out patients with no fio2 measurements (3)

length(unique(amd$admissionid)) #1408

fio2_measured <- amd %>% filter(var == "fio2") %>%
  distinct(admissionid) # 1406 so 2 additional patient had no fio2 measurement after filtering out observations concurrent with IMV

temp <- unique(amd$admissionid)
temp[!(temp %in% unique(fio2_measured$admissionid))]
# 8830 13788 have no fio2 measurements after filtering

amd <- filter(amd,
              admissionid %in% fio2_measured$admissionid)

# exclude some patients with first fio2 < 40 that snuck in due to validated / unvalidated
fio2geq40 <- 
  amd %>%
  group_by(admissionid) %>%
  filter(chartminute <= 0) %>%
  filter(var == "fio2") %>%
  summarise(fio2 = last(value, order_by = chartminute)) %>%
      filter(fio2 >= 40) # 1 patient, 7708

temp <- unique(amd$admissionid)
temp[!(temp %in% unique(fio2geq40$admissionid))]


amd <- filter(amd,
              admissionid %in% fio2_measured$admissionid) %>%
  filter(admissionid %in% fio2geq40$admissionid)


# convert pressor measurements to binary on / off

pressor <- 
  amd %>%
    filter(grepl("Noradrenaline", var) |
             grepl("Adrenaline", var) |
             grepl("Dopamine", var)) %>%
    group_by(admissionid, var, chartminute) %>%
    summarise(pressor = sum(value)) %>%
    mutate(pressor = ifelse(pressor == 0, NA, pressor)) %>%
    mutate(value = ifelse(pressor < 0, 0, pressor),
           var = "pressor") %>%
    group_by(admissionid) %>%
    fill(value) %>%
  select(-pressor) %>%
  distinct()

amd <- filter(amd, var %in% c("fio2",
                              "GCS",
                              "heart_rate",
                              "resp_rate",
                              "mbp",
                              "o2_device",
                              "pco2",
                              "po2",
                              "ph",
                              "spo2")) %>%
  bind_rows(pressor)



# remove tracheostomy patients
# 4, 18, 19, 11, 13, 14, 15 are all BVM or tracheostomy related

# o2 device counts
amd %>% filter(var == "o2_device") %>%
  mutate(value =  factor(value)) %>%
  group_by(admissionid, value, .drop = F) %>%
  summarise(measured = n() > 0) %>%
  group_by(value) %>%
  summarise(out = sum(measured))

trach <- amd %>% filter(var == "o2_device") %>%
  mutate(value =  factor(value)) %>%
  group_by(admissionid, value, .drop = F) %>%
  summarise(measured = n() > 0) %>%
  group_by(admissionid) %>%
  filter(value %in% c(4, 18:19, 11, 13:15)) %>%
  summarise(out = sum(measured==T)>0) %>%
  filter(out == T)

nrow(trach) # 25
# filter them out

amd <- amd %>%
  filter(!(admissionid %in% trach$admissionid))

imv_ever <- imv_ever %>%
  filter(!(admissionid %in% trach$admissionid)) %>%
  filter(admissionid %in% fio2_measured$admissionid)

patients <- amd %>%
  group_by(admissionid) %>%
  summarise(admissionid = first(admissionid)) %>%
  left_join(imv_ever, by = "admissionid") %>%
  # filter out trach patients
  filter(!(admissionid %in% trach$admissionid)) %>%
  # join with baseline data
  left_join(amd_b, by = "admissionid") %>%
  select(patientid, admissionid, location, admissionyeargroup,
         admissioncount,
         dischargedat, gender, agegroup, dateofdeath, weightgroup, heightgroup, 
         specialty, imv_time) 
  

# make sure only one admission per patient (take the first one)
patients %>%
  group_by(patientid) %>%
  summarise(count = n()) %>%
  filter(count > 1) %>%
  ungroup() %>% 
  summarise(n()) # 84 admits from 78 pts

patients <- patients %>%
  group_by(patientid) %>%
  mutate(firstadmit = min(admissioncount)) %>% 
  filter(admissioncount == firstadmit) %>%
  select(-admissioncount, -firstadmit)

amd <- filter(amd, 
              admissionid %in% patients$admissionid)

# and one patient has only negative times so remove them
filter(amd, admissionid == 18674)

patients <- filter(patients, !(admissionid == 18674))
amd <- filter(amd, !(admissionid == 18674))

# break o2 devices into categories
# 16 = NIV
# 17 = NRB
# 20 = IMV
# all others coded as all 0s for the above

niv <- amd %>%
  filter(var == "o2_device") %>%
  mutate(var = "niv",
         value = ifelse(value == 16, 1, 0)) %>%
  distinct()

stdo2 <- amd %>%
  filter(var == "o2_device") %>%
  mutate(var = "stdo2",
         value = ifelse(value == 17, 1, 0)) %>%
  distinct()

amd <- amd %>%
  filter(!(var == "o2_device")) %>%
  bind_rows(niv) %>%
  bind_rows(stdo2)


amd %>%
  group_by(var) %>%
  summarise(min = min(value, na.rm = T),
            max = max(value, na.rm = T))

  

# some duplicate measurements at the same time for the same patient

amd %>% 
  group_by(admissionid, var, chartminute) %>%
  filter(!is.na(value)) %>%
  count() %>%
  filter(n > 1)
# 876 / 348849 = 0.25%

cont_vars <- c("heart_rate",
               "mbp",
               "po2",
               "pco2",
               "ph",
               "resp_rate",
               "spo2",
               "fio2",
               "GCS")

# take the mean (continuous) and max (binary) because no other basis to select between them.
# fortunately a small proportion
amd <- amd %>% 
  group_by(admissionid, var, chartminute) %>%
  filter(!is.na(value)) %>%
  summarise(value = ifelse(var %in% cont_vars,
                           mean(value),
                           max(value))) %>%
  ungroup() %>%
  select(admissionid, chartminute, var, value) %>%
  arrange(admissionid, chartminute)

# transform, center, scale

log_group <- c("heart_rate",
               "mbp",
               "po2",
               "pco2",
               "ph",
               "resp_rate")

transform <- function(row, log_group){
  var <- row[3]
  value <- as.numeric(row[4])
  if(var %in% log_group){log(value)}
  # use logit transformation functions for GCS, spo2, and fio2
  # because they have a restricted range 
  else if (var == "GCS"){logit((value-2.9)/12.2)}
  else if (var == "spo2"){logit(value/101)} # take the floor when transforming back
  else if (var == "fio2"){logit((value-20.9)/79.3)} # take the floor when transforming back
  else value
}

amd_transformed <- amd
amd_transformed$value <- apply(amd, 1, transform, log_group)

cont_vars <- c("heart_rate",
               "mbp",
               "po2",
               "pco2",
               "ph",
               "resp_rate",
               "spo2",
               "fio2",
               "GCS")

amd_tcs <- amd_transformed

mean_sd <- filter(amd_tcs, var %in% cont_vars) %>%
  group_by(var) %>%
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T))

saveRDS(mean_sd, file = "amd_oc_mean_sd.rds")

for (j in 1:nrow(mean_sd)){
  var2 <- as.character(mean_sd[j,1])
  mean <- as.double(mean_sd[j,2])
  sd <- as.double(mean_sd[j,3])
  
  amd_tcs <- mutate(
    amd_tcs,
    value = ifelse(var == var2,
                   (value-mean)/sd,
                   value))
}


# obs durations per person
lastobs <- amd %>%
  group_by(admissionid) %>%
  summarise(lastobs = max(chartminute))

# prep for stan model and interpolation

amd_tcs_stan <- filter(amd_tcs,
                       var %in% c("resp_rate",
                                       "spo2",
                                       "heart_rate",
                                      # "mbp",
                                       "fio2",
                                       "po2",
                                       "ph",
                                       "pco2",
                                       "GCS",
                                       "pressor",
                                       "niv",
                                       "stdo2"),
                       chartminute > -12*60) %>% 
  rename(time = chartminute) %>%
  distinct() %>%
  filter(!is.na(value)) %>%
  pivot_wider(id_cols = c(admissionid,time),
              names_from = var,
              values_from = value,
              values_fill = NA) %>%
  group_by(admissionid) %>%
  fill(fio2, pressor, niv, stdo2) %>%
  pivot_longer(cols = heart_rate:GCS,
               names_to = "var",
               values_to = "value") %>%
  filter(time >= 0 | !(is.na(value))) %>%
  mutate(value = ifelse(var == "pressor" & is.na(value),
                        0, value)) %>%
  filter(!(var %in% c("ph","pco2","po2","fio2","niv","stdo2") &
             is.na(value)))

amd_tcs_stan %>%
  group_by(var) %>%
  summarise(nasum = sum(is.na(value)),
            namean = mean(is.na(value)),
            count = sum(!is.na(value)))%>%
  write_csv(file = "amd_missing.csv")

amd_tcs_stan %>%
  filter(!is.na(value)) %>%
  group_by(admissionid, var) %>%
  summarise(count = n()) %>%
  left_join(lastobs, by = "admissionid") %>%
  mutate(lastobs = lastobs+1) %>%
  group_by(var) %>%
  summarise(mean(lastobs)/mean(count)) %>%
  write_csv(file = "amd_time_btw_obs.csv")


#########################################

# Split cohort into manageable bites
# to find Gaussian process hyperparameters
# take a random contiguous 48h from every patient

amd_b <- patients

K_total <- nrow(amd_b)
set.seed(20210205)
interval_length <- 48*60
K_subset <- sample(1:K_total, size = 400, replace = F)
patients <- data.frame(id = 1:K_total, 
                       admissionid = amd_b$admissionid)

patients2 <- filter(patients, id %in% K_subset) %>%
  mutate(id = 1:length(K_subset))

amd_b <- left_join(amd_b, patients, by = "admissionid")

library(cmdstanr)

# prep data for stan model

# data for interpolating each patient

cvars <- c("resp_rate",  "heart_rate", "spo2","fio2", 
           "GCS", "po2", "ph","pco2") # continuous variables 
bvars <- c("pressor", "niv","stdo2")

D_c <- length(cvars)
D_b <- length(bvars)
D = D_c + D_b

dft <- left_join(patients, amd_tcs_stan, by = "admissionid") %>%
  mutate(var = factor(var, levels = c(cvars,bvars), labels = 1:D)) %>%
  mutate(time = time/(24*60)) %>%
  arrange(id, var, time) %>%
  mutate(var = as.integer(as.character(var))) 

dft_pred <- filter(dft, is.na(value)) %>%
  filter(var %in% c(1:3,5)) %>%
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

saveRDS(interp_dat, file = "C:/Users/chris/OneDrive/Documents/LargeDataFiles/amd_interpolation_dat.rds")

amd_b <- 
  amd_b %>%
    left_join(amd_elig, by = "admissionid") %>%
    mutate(dateofdeath = dateofdeath/60000-eligibletime)
    

saveRDS(amd_b, "oc_amd_b.rds")
saveRDS(amd, "oc_amd.rds")
saveRDS(amd_tcs, "oc_amd_tcs.rds")
saveRDS(amd_tcs_stan, "oc_amd_tcs_stan.rds")

# version to fit initial GP

dft <- left_join(patients2, select(dft,-id), by = "admissionid") %>%
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

write_stan_json(data = interp_GP_dat, file = "C:/Users/chris/OneDrive/Documents/LargeDataFiles/amd_data_for_GP_interpolation.json")


# fixed forward fill imputation

order <- data.frame(var = c("resp_rate",  "heart_rate", "spo2","fio2",  "GCS",
                            "po2","ph","pco2",
                            "pressor", "niv", "stdo2"),
                    num = 1:D)

vars <- left_join(order, mean_sd, by = "var")


amd_fi <- 
  amd_tcs_stan %>%  
  pivot_wider(names_from = "var",
              values_from = "value") %>%
  fill(c(heart_rate:spo2, resp_rate:GCS)) %>%
    mutate(pressor = ifelse(is.na(pressor),0,pressor)) %>%
  pivot_longer(cols = heart_rate:GCS,
               names_to = "var") %>%
  filter(time >= 0, !is.na(value)) %>%
  left_join(vars, by = "var") %>%
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
  mutate(time = time/60) %>%
  select(-mean, -sd)

saveRDS(amd_fi, file = "oc_amd_fi.rds")
