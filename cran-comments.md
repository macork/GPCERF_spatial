Re-submission (January 21, 2023)

Thank you for taking the time to review the GPCERF 0.2.0 package. I fixed the 
URL issues.  


Best regards,
Naeem Khoshnevis
RCD - HUIT Harvard University



Submission (January 21, 2023)
Thank you for taking the time to review the GPCERF 0.2.0 package. In this update we:
- Enhanced both approaches to use multiple CPU cores.
- Updated default values when necessary and removed default values from some other functions. 
- Converted a couple of functions into internal functions (the package does not have reverse dependencies.)
- Reviewed the code to follow the tidyverse styling guide.
- Improved test coverage.

We also tested the package on numerous builds, including:

- Windows Server 2022, R-release, 32/64 bit (windows-x86_64-release)
- macOS 10.13.6 High Sierra, R-release, brew (macos-highsierra-release)
- Debian Linux, R-devel, clang, ISO-8859-15 locale (debian-clang-devel)
- Ubuntu Linux 20.04.1 LTS, R-devel, GCC (ubuntu-gcc-devel)
- Fedora Linux, R-devel, clang, gfortran (fedora-clang-devel)


Best regards,
Naeem Khoshnevis
RCD - HUIT Harvard University
