# pcm_code
Code for the simulations in the paper "The Projected Covariance Measure for assumption-lean variable significance testing".

## Overview
The package contains the following files:

#### `test_functions.R`
This file contains all the functions used to compute the PCM test and the tests that we compare to.

#### `sim_gam_binary_comparison.R`


#### `sim_gam_comparison.R`


#### `sim_linear_rates.R`


#### `sim_ranger_comparison.R`


#### `sim_regression_errors.R`


## Applying test functions with other regression methods




## Package versions
Below is included the `sessionInfo()` output of an R session where the code is known to be working. It is most crucial to ensure that the `vimp` package has the correct version as the interface of the `cv_vim` function has changed considerably over time. 

```{r}
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 13.5

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] mgcv_1.8-41   nlme_3.1-160  vimp_2.3.1    ranger_0.14.1

loaded via a namespace (and not attached):
 [1] shape_1.4.6           tidyselect_1.2.0      dread_0.0.0.9000      kernlab_0.9-31        reshape2_1.4.4        splines_4.2.0         lattice_0.20-45      
 [8] RPtests_0.1.5         vctrs_0.6.3           generics_0.1.3        lawstat_3.5           utf8_1.2.3            survival_3.4-0        rlang_1.1.1          
[15] pracma_2.4.2          pillar_1.9.0          glue_1.6.2            xgboost_1.6.0.1       RColorBrewer_1.1-3    quantregForest_1.3-7  foreach_1.5.2        
[22] lifecycle_1.0.3       mize_0.2.4            plyr_1.8.8            stringr_1.5.0         grf_2.2.1             caTools_1.18.2        mvtnorm_1.1-3        
[29] codetools_0.2-18      fansi_1.0.4           Rcpp_1.0.9            jsonlite_1.8.3        Kendall_2.2.1         stringi_1.7.8         gam_1.22-2           
[36] dplyr_1.1.2           CondIndTests_0.1.5    rbibutils_2.2.10      grid_4.2.0            Rdpack_2.4            cli_3.6.1             tools_4.2.0          
[43] bitops_1.0-7          magrittr_2.0.3        glmnet_4.1-4          tibble_3.2.1          randomForest_4.7-1.1  weightedGCM_0.1.0     pkgconfig_2.0.3      
[50] MASS_7.3-58.1         Matrix_1.5-3          data.table_1.14.8     SuperLearner_2.0-28.1 nnls_1.4              rstudioapi_0.14       iterators_1.0.14     
[57] R6_2.5.1              boot_1.3-28.1         compiler_4.2.0 
```
