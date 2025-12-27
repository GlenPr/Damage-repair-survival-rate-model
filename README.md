# Damage-repair-survival-rate-model
For modelling multivariate longitudinal binary data with survival, for example health deficits.

# Purpose
To analyze transitions and their dependence on age and the frailty index (FI; defined as average of binary values). Original purpose was to look for changes in robustness and resilience. It ultimately found a tipping point.

# Background
The FI is defined as the average over all binary values an individual has (at a particular time point). Damage is defined by the transition 0->1. Repair is defined by the transition 1->0. Survival is a first passage event e.g. death.

Robustness is defined as the ability to resist damage transitions while resistance is the ability to repair. Research in the field of geroscience has found evidence that both robustness and resilience decrease with age and increasing frailty index. We sought to check this quantitatively using this package.

# Installation
Download and open in R (no install necessary). The functions are ready to go (.R). The vignette works in RStudio (.Rmd). Works in R v4.2.2 and RStudio 2023.03.0. Tested on Windows 11.

# Dependencies
-survival

-ggplot2

-metr (recommended)

-ggridges (recommended)

-mgcv (recommended)

# What it does
Takes multivariate binary longitudinal data. Originally used to analyse binary health deficits (0: good, 1: bad). Estimates parameters for a log-rate model where the damage, repair and survival rates depend on the current FI and age of the individual. The model must be linear in age but can have arbitrary dependence on the FI; use GenBasis() to determine the number and form of the FI dependence.

# How do I use it?
Fitting is entirely self-contained in a single function that will bootstrap and generate parameter estimates, errors and diagnostic plots. Please start with the vignette, it explains dependenceis, how to shape data, and perform the analysis. Runtime is about 5 minutes total. 

# For more information / cite as
Pridham, G., Rockwood, K. & Rutenberg, A. D. Dynamical modelling of the frailty index indicates that health reaches a tipping point near age 75. arXiv [q-bio.QM] (2024).
  

@ARTICLE{Pridham2024-su,
  title         = "Dynamical modelling of the frailty index indicates that
                   health reaches a tipping point near age 75",
  author        = "Pridham, Glen and Rockwood, Kenneth and Rutenberg, Andrew D",
  journal       = "arXiv [q-bio.QM]",
  month         =  dec,
  year          =  2024,
  url           = "https://scholar.google.com/citations?view_op=view_citation&hl=en&citation_for_view=cVHiV0gAAAAJ:ZeXyd9-uunAC",
  archivePrefix = "arXiv",
  primaryClass  = "q-bio.QM",
  eprint        = "2412.07795"
}

