Resubmission (June 30, 2022)

Thank you so much for taking the time and reviewing the GPCERF 0.1.0 package and 
providing feedback.

Here are the changes based on your recommendations:

- Converted the reference in the DESCRIPTION file into authors (year)<arXiv:...>
format, and put the title in the quotes. 

- Removed examples from the internal (not-exported) functions:
  - calc_ac
  - compute_deriv_nn
  - compute_sd_gp

- The example in the estimate_cref_nngp function takes longer than 5s to execute.
Thanks for pointing that out. I replaced \dontrun{} with \donttest{}.



Best regards, 
Naeem Khoshnevis 
FASRC - Harvard University


Resubmission (June 28, 2022)

Thank you so much for taking the time and reviewing the GPCERF 0.1.0 package. 

I fixed the url and updated the site. 

Best regards, 
Naeem Khoshnevis 
FASRC - Harvard University



Original Submission (June 28, 2022)

Thank you so much for taking the time and reviewing the GPCERF 0.1.0 package. 

This R package, Provides a non-parametric Bayesian framework based on Gaussian process priors for estimating the causal effects of continuous exposure and detecting change points in the causal exposure response curves using observational data.

Here is the method paper which is under review:

Ren, B., Wu, X., Braun, D., Pillai, N. and Dominici, F., 2021. Bayesian 
modeling for exposure response curve via Gaussian processes: Causal effects of 
exposure to air pollution on health outcomes. arXiv preprint arXiv:2105.03454.

We followed the best practices in the package development and successfully tested the package on the following rhub platforms:

1) macOS 10.13.6 High Sierra, R-release, brew (macos-highsierra-release)
2) Debian Linux, R-devel, GCC ASAN/UBSAN (linux-x86_64-rocker-gcc-san)
3) Debian Linux, R-devel, clang, ISO-8859-15 locale (debian-clang-devel)
4) Ubuntu Linux 20.04.1 LTS, R-devel, GCC (ubuntu-gcc-devel)
5) Windows Server 2022, R-devel, 64 bit (windows-x86_64-devel)

Sometimes, I get the following note: 

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
  
There is a discussion in the following Github issue about this. 

https://github.com/r-hub/rhub/issues/503 

It seems like a bug with miktex. Please let me know if you have any suggestions for fixing it. 


Best regards, 
Naeem Khoshnevis 
FASRC - Harvard University
