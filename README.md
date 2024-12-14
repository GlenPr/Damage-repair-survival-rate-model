# Damage-repair-survival-rate-model
For modelling multivariate longitudinal binary data with survival, for example health deficits.

# Purpose
To analyze transitions and their dependence on age and the frailty index (FI; defined as average of binary values). Original purpose was to look for tipping point behaviour.

# Background
The FI is defined as the average over all binary values an individual has (at a particular time point). Damage is defined by the transition 0->1. Repair is defined by the transition 1->0. Survival is a first passage event e.g. death.

# Installation
Download and open in R. The functions are ready to go (.R). The vignette works in R studio (.Rmd).

# Dependencies
-survival

-ggplot2

-metr (recommended)

-ggridges (recommended)

-mgcv (recommended)

# What it does
Takes multivariate binary longitudinal data. Originally used to analyse binary health deficits (0: good, 1: bad). Estimates parameters for a log-rate model where the damage, repair and survival rates depend on the current FI and age of the individual. The model must be linear in age but can have arbitrary dependence on the FI; use GenBasis() to determine the number and form of the FI dependence.

# How do I use it?
Fitting is entirely self-contained in a single function that will bootstrap and generate parameter estimates, errors and diagnostic plots.

# For more information / cite as
Pridham G, Rockwood K, Rutenberg AD. Dynamical modelling of the frailty index indicates that health reaches a tipping point near age 75. arXiv preprint arXiv:2412.07795. 2024 Dec 2.

@article{pridham2024dynamical,
  title={Dynamical modelling of the frailty index indicates that health reaches a tipping point near age 75},
  author={Pridham, Glen and Rockwood, Kenneth and Rutenberg, Andrew D},
  journal={arXiv preprint arXiv:2412.07795},
  year={2024}
}
