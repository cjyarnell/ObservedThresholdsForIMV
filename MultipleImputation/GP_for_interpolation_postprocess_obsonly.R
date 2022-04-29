
library(cmdstanr)
library(dplyr)
library(bayesplot)
library(ggpubr)
library(posterior)
library(reshape2)
library(ggplot2)
library(data.table)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

options(mc.cores=6)

# colours
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_blue = "cornflower blue"

# inv logit
inv_logit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

D <- 9

order <- data.frame(variable = c("resp_rate",  "heart_rate", "spo2","fio2",  
                                 "GCS", "wob", "hfnc", "niv", "stdo2"),
                    var = 1:D)

# visualize the covariance matrix
visualize_cov <- function(X, n_iter, n_chain){
  X <- array(as.numeric(unlist(X)), dim = c(n_iter, n_chain, D, D))
  cov <- array(dim = c(n_iter*n_chain, D, D))
  pos = 1
  for (n in 1:n_iter){
    for (m in 1:n_chain){
      cov[pos,,] <- tcrossprod(X[n,m,,])
      #      colnames(cov[pos,,]) <- vars$variable
      #      rownames(cov[pos,,]) <- vars$variable
      pos = pos + 1
    }
  }
  cov_mean <- apply(cov, c(2,3), mean)
  cov_mean[cov_mean>1] <- 1
  
  melted_cormat <-reshape2::melt(cov_mean)
  
  mc1 <- melted_cormat %>%
    left_join(order, by = c("Var1" = "var")) %>%
    left_join(order, by = c("Var2" = "var")) %>%
    mutate(value = ifelse(Var1 < Var2,
                          0, value)) %>%
    mutate(Var1 = factor(Var1,
                         labels = order$variable),
           Var2 = factor(Var2,
                         labels = order$variable))
  
  
  mc2 <- melted_cormat %>%
    left_join(order, by = c("Var1" = "var")) %>%
    left_join(order, by = c("Var2" = "var")) %>%
    mutate(value = ifelse(Var1 >= Var2,
                          0, value)) %>%
    mutate(value = paste0(round(value, 2))) %>%
    mutate(value = ifelse(Var1 >= Var2,
                          "", value),
           Var1 = factor(Var1,
                         labels = order$variable),
           Var2 = factor(Var2,
                         labels = order$variable)) 
  
  # Heatmap
  plot <- ggplot(data = mc1, aes(y = Var2,
                                 x = Var1))+
    geom_tile(color = "white",
              aes(fill = value))+
    geom_text(data = mc2,
              aes(label = value),
              size = 2,
              color = "grey50") + 
    scale_fill_gradient2(low = c_dark, high = c_blue, mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Correlation",
                         breaks = c(-1,0,1)) +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 90))+
    labs(y = "", x = "", title = "Covariance between time series variables") +
    coord_fixed()
  
  list(cov, plot)
}


# MCMC fit

files = c("models/interpolation_obsonly_GP1.csv",
         "models/interpolation_obsonly_GP2.csv",
         "models/interpolation_obsonly_GP3.csv",
         "models/interpolation_obsonly_GP4.csv")

interpolation_GPfit <- read_cmdstan_csv(files,
                                        variables = c("rho",
                                                "alpha",
                                                "sigma",
                                                "offset",
                                                "L_C"))


# trace plots
rho_trace <- mcmc_trace(interpolation_GPfit$post_warmup_draws, pars = "rho")
alpha_trace <- mcmc_trace(interpolation_GPfit$post_warmup_draws, regex_pars = "alpha")
sigma_trace<- mcmc_trace(interpolation_GPfit$post_warmup_draws, regex_pars = "sigma")
offset_trace<- mcmc_trace(interpolation_GPfit$post_warmup_draws, regex_pars = "offset")

svg("plots/interpolation_GP_rho_trace_obsonly.svg",
   width = 5, height = 5)
rho_trace
dev.off()

svg("plots/interpolation_GP_alpha_trace_obsonly.svg",
   width = 10, height = 8)
alpha_trace
dev.off()

svg("plots/interpolation_GP_sigma_trace_obsonly.svg",
   width = 10, height = 5)
sigma_trace
dev.off()

svg("plots/interpolation_GP_offset_trace_obsonly.svg",
   width = 5, height = 5)
offset_trace
dev.off()

fit <- interpolation_GPfit

meanrho <- mean(fit$post_warmup_draws[,,grepl("rho", dimnames(fit$post_warmup_draws)$variable)])
meanalpha <- apply(fit$post_warmup_draws[,,grepl("alpha", dimnames(fit$post_warmup_draws)$variable)], 
                   3, mean)
meansigma <- apply(fit$post_warmup_draws[,,grepl("sigma", dimnames(fit$post_warmup_draws)$variable)], 
                   3, mean)
meanoffset <- apply(fit$post_warmup_draws[,,grepl("offset", dimnames(fit$post_warmup_draws)$variable)], 
                   3, mean)
meanL_C <- matrix(apply(fit$post_warmup_draws[,,grepl("L_C", dimnames(fit$post_warmup_draws)$variable)], 
                        3, mean), nrow = D, ncol = D)

n_chain = 4
n_iter = 200
Cov_C <- visualize_cov(fit$post_warmup_draws[,,grepl("L_C", dimnames(fit$post_warmup_draws)$variable)], 
                       n_iter, n_chain)

svg("plots/interpolation_GP_cov.svg", width = 7, height = 7)
Cov_C[[2]]
dev.off()

GP_params <- list(rho = meanrho,
                  alpha = meanalpha,
                  sigma = meansigma,
                  offset = meanoffset,
                  L_C = meanL_C)
saveRDS(GP_params, "data/gp_params_obsonly.rds")
