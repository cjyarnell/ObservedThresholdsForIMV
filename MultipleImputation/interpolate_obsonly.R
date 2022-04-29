#!/usr/bin/Rscript
# this program interpolates one patient's trajectories
# and saves them as an .rds file

setwd("/gpfs/fs0/scratch/g/gtomlins/cyarnell/ObservedCriteria")
library(dplyr)

#start <- as.integer(commandArgs(trailingOnly = TRUE))
#cycle <- 897

# colours
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_blue = "cornflower blue"

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

library(cmdstanr)
library(tidyr)
library(ggplot2)
library(dplyr)
source("/gpfs/fs0/scratch/g/gtomlins/cyarnell/ObservedCriteria/models/check_thresholds_obsonly.R")


# load the interpolation data

dat <- readRDS("./data/interpolation_dat_obsonly.rds")

# completed patients
completed <- list.files("interpolations2/")
completed <- data.frame(filename = completed) %>% 
  mutate(stay_id = as.numeric(gsub(".rds","", gsub("interpolations_","",filename))))

patients <- dat$patients %>%
    left_join(completed, by = "stay_id") %>%
    filter(is.na(filename))

GP_params <- readRDS("data/gp_params_obsonly.rds")


interpolate <- function(pt_stay_id){
    
print(paste0("starting interpolation for ", pt_stay_id))

id = filter(patients, stay_id == pt_stay_id)$id
    
# filter to single patient
k_obs <- dat$k_obs

N_obs <- sum(k_obs == id)
D <- dat$D
D_c <- dat$D_c
D_b <- dat$D_b
x_obs <- dat$x_obs[k_obs == id]
d_obs <- dat$d_obs[k_obs == id]
d_obs_matrix <- dat$d_obs_matrix[k_obs == id,]

N_c <- dat$N_kc[id]
N_b <- dat$N_kb[id]

start_c <- ifelse(id == 1, 0, sum(dat$N_kc[1:(id-1)]))
start_b <- ifelse(id == 1, 0, sum(dat$N_kb[1:(id-1)]))

y_c <- array(dat$y_c[(start_c+1):(start_c + N_c)])
y_b <- array(dat$y_b[(start_b+1):(start_b + N_b)])

k_pred <- dat$k_pred
x_pred <- array(dat$x_pred[k_pred == id])
d_pred <- array(dat$d_pred[k_pred == id])
N_pred <- sum(k_pred == id)

temp_dat <- list(N_obs = N_obs,
                 D = D,
                 x_obs = x_obs,
                 d_obs = d_obs,
                 d_obs_matrix = d_obs_matrix,
                 
                 N_c = N_c,
                 D_c = D_c,
                 y_c = y_c,
                 
                 N_b = N_b,
                 D_b = D_b,
                 y_b = y_b,
                 
                 N_pred = N_pred,
                 x_pred = x_pred,
                 d_pred = d_pred
                )

# combine with fixed GP parameters


temp_dat <- c(temp_dat, GP_params)

# run the stan fit

warmup = 200
sample = 100

interp_mod <- cmdstan_model("models/multi_interp_fixedGP_obsonly.stan")

interp_fgp <- interp_mod$sample(data = temp_dat,
                                iter_warmup = warmup,
                                iter_sampling = sample,
                                chains = 1,
                               refresh = 5)

# extract the predictions

chains <- 1

y_pred <- interp_fgp$draws(variables = "y_pred")[,1,]

df_pred <- data.frame(stay_id = pt_stay_id,
                      time = rep(x_pred, each = sample),
                      iter = rep(1:sample, times = N_pred),
                      var = rep(d_pred, each = sample),
                      value = y_pred[1:(N_pred*sample)])

# combine with observed data
    
df_obs <- data.frame(stay_id = pt_stay_id,
                     time = rep(x_obs, each = sample),
                     iter = rep(1:sample, times = N_obs),
                     var = rep(d_obs, each = sample),
                     value = rep(c(y_c, y_b), each = sample))
    
df <- bind_rows(df_obs, df_pred) %>%
    arrange(time, var) %>%
    # convert time back to hours
    mutate(time = time*24)
    
# reverse the transformations

mean_sd <- readRDS("data/oc_mean_sd2.rds") %>%
  filter(variable %in% c("resp_rate",
                         "spo2",
                         "GCS",
                         "fio2",
                         "heart_rate",
                         "po2",
                         "ph",
                         "pco2"))

order <- data.frame(variable = c("sbp", "resp_rate",  "heart_rate", "spo2","fio2",  "GCS",
                                 "po2","ph","pco2",
                                 "wob", "pressor", "hfnc", "niv", "stdo2"),
                    var = 1:D)

vars <- left_join(order, mean_sd, by = "variable")

df <- left_join(df, vars, by = "var") %>%
  mutate(value = ifelse(var <= D_c, value*sd + mean, value))

log_group <- c("resp_rate",
               "heart_rate",
              "sbp",
              "po2",
              "ph",
              "pco2")

df <- mutate(df,
                  value = ifelse(variable %in% log_group,
                                 exp(value),
                                 value)) %>% 
  mutate(value = ifelse(variable == "spo2",
                        round((inv_logit(value))*101),
                        value)) %>%
  mutate(value = ifelse(variable == "GCS",
                        round((inv_logit(value))*12.2+2.9),
                        value)) %>%
  mutate(value = ifelse(variable == "fio2",
                        round((inv_logit(value))*79.3+20.9),
                        value))

# get the baseline characteristics

pt_bl <- readRDS("data/oc_pdx_b.rds") %>%
  filter(stay_id == pt_stay_id)

df <- filter(df, time >= 0)
    
interpolation_list <- list(stay_id = pt_stay_id,
                           pt_bl = pt_bl,
                           df = df)

saveRDS(interpolation_list,
        file = paste0("interpolations2/interpolations_",
                      pt_stay_id,
                      ".rds"))

print(paste0("finished patient ", id))
}

library(doParallel)

print("start parallelization")

ncores = Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=80)

output <- foreach(i = 1:(nrow(patients))) %dopar% interpolate(patients$stay_id[i])