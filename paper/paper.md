---
title: 'GPCERF - An R package for implementing Gaussian processes for estimating causal exposure response curves'
tags:
  - R
  - causal inference 
  - Gaussian Processes
  - causal exposure response function
authors:
  - name: Naeem Khoshnevis^[Corresponding author]
    orcid: 0000-0003-4315-1426
    affiliation: "1"
  - name: Boyu Ren
    orcid: 0000-0002-5300-1184
    affiliation: "2"
  - name: Danielle Braun
    orcid: 0000-0002-5177-8598
    affiliation: "3"
affiliations:
 - name: Research Computing, Harvard University, Cambridge, MA, United States of America
   index: 1
 - name: McLean Hospital, Belmont, MA, United States of America
   index: 2
 - name: Department of Biostatistics, Harvard School of Public Health, Cambridge, MA, United States of America
   index: 3

date: 15 March 2023
bibliography: refs.bib
---

# Summary

We present the GPCERF R [@R_2022] package, which employs a novel Bayesian approach based on Gaussian Process (GP) to estimate the causal exposure-response function (CERF) for continuous exposures, along with associated uncertainties. R packages that target at these causal effects under a binary exposure setting exist @[citep][MatchIt_R], as well as in the continuous exposure setting @[citep][CausalGPS_R]. However, they often rely on a separate resampling stage to quantify uncertainty of the estimates, which can be computationally demanding. GPCERF provides a two-step end-to-end solution for causal inference with continuous exposures that is equipped with automatic and efficient uncertainty quantification. During the first step (design phase), the algorithm searches for optimal hyperparameters (using the exposures and covariates) that yield the optimal covariate balance, allowing for the creation of a balanced data set. The selected hyperparameters are then used in the second step (analysis phase) to estimate the CERF on the balanced data set and its associated uncertainty using two different types of GPs: a standard GP and a nearest-neighbor GP (nnGP). Standard GP offers high accuracy in estimating CERF but is also computationally intensive. The nnGP is a computationally efficient approximation of the standard GP and is well-suited for the analysis of large-scale dataset. 
