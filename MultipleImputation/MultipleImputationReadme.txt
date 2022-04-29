This folder contains the programs used to:

1) Generate the Gaussian process model for multiple imputation (MIMIC: GP_for_interpolation.stan, AmsterdamUMCdb: amds_GP_for_interpolation.stan)
2) Extract the mean hyperparameters from the Gaussian process models (_postprocess.R files)
3) Run the scripts that impute missing data (interpolate_obsonly.R and amd_interpolate.R) which are parallelized by patient