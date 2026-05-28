## Resubmission
This is a resubmission. In this version I have:

* Added a reference to the source paper in the DESCRIPTION file.

* Removed print functions from sample_partition_correlation() and sample_partition_independence() functions 
  and replaced them with stop functions that include more informative error messages.
  
* Changed "suppressPrint=FALSE" option in mcmc_raceNMA() function to "verbose=TRUE", to better comply with R style.

* Updated reproducibility vignette accordingly.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
