This folder contains the code used to calculate the outcomes (ie to see when thresholds are met and what happens after those thresholds are met)

The "check_thresholds" files contain the functions for checking whether or not a threshold has been met.
The "outcomes" files calculate whether or not a threshold was met and if so, whether or not invasive ventilation followed, for all patients across all imputations.
The "obscrit_categorical.stan" file contains the code for the adjusted analysis, with categorical outcomes whose latent outcome variables are correlated by a linear covariance structure.
