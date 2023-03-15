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

We present the GPCERF R [e.g., @R_2022] package, which employs a novel Bayesian approach based on Gaussian Process (GP) to estimate the causal exposure-response function (CERF) for continuous exposures, along with associated uncertainties. R packages that target at these causal effects under a binary exposure setting exist [e.g., @MatchIt_R], as well as in the continuous exposure setting [@CausalGPS_R]. However, they often rely on a separate resampling stage to quantify uncertainty of the estimates, which can be computationally demanding. GPCERF provides a two-step end-to-end solution for causal inference with continuous exposures that is equipped with automatic and efficient uncertainty quantification. During the first step (design phase), the algorithm searches for optimal hyperparameters (using the exposures and covariates) that yield the optimal covariate balance, allowing for the creation of a balanced data set. The selected hyperparameters are then used in the second step (analysis phase) to estimate the CERF on the balanced data set and its associated uncertainty using two different types of GPs: a standard GP and a nearest-neighbor GP (nnGP). Standard GP offers high accuracy in estimating CERF but is also computationally intensive. The nnGP is a computationally efficient approximation of the standard GP and is well-suited for the analysis of large-scale dataset. 

# Statement of need

Existing R packages for estimating causal exposure response curve of continuous exposures usually requires resampling approaches, such as boostrap, to obtain uncertainty of the estimates [e.g., @CausalGPS_R]. However, when the number of observations is large, resampling based algorithms can be computationally expensive. To address this gap, we have implemented a novel Bayesian approach that uses Gaussian Processes (GPs) prior for the counterfactual outcome surfaces to enable flexible estimation of CERF. By leveraging the posterior distributions of the model parameters, we can automatically quantify the uncertainty of the estimated CERF [for more details see @Ren_2021_bayesian]. Since standard GPs are notorious for its lack of scalability due to the operations involving inversion of the covariance matrices, we have also added in our approach a nearest-neighbour GP (nnGP) prior to achieve computationally efficient inference of CERF when dealing with large-scale datasets.

# Overview

In the context of causal inference for continuous exposure, the main target for inference is the so-called casual exposure response function (CERF), which is defined as the expectation of the counterfactual outcomes over a range of exposure levels in a given population. If we denote the counterfactual outcome at exposure level $w$ by $Y(w)$, CERF is indeed the function $R(w) = \mathbb E(Y(w))$ defined on a set of $w$ of interest. One should be careful when extracting CERF from observational data, which usually contain not only the outcome $Y$ and exposure $W$, but also potential confounders $C$. The main reason is that the confounders $C$ leads to biased inference of $R(w)$ if not adjusted properly, such as through a simple regression model. We choose to follow the approach in @Hirano_2004, which is based on generalized propensity score (GPS), in this package to handle this issue. GPS, denoted by $s(W,C)$, is defined as the conditional density of $W$ given $C$. It has been shown that one can obtain an unbiased estimator of the causal effect of $W$ provided the conditional distribution of $Y$ given $W$ and $s(W,C)$ is known @Hirano_2004.

# References
