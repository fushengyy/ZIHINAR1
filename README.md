
<!--  README.md is generated from README.Rmd. Please edit that file -->

The First-order Integer-valued Autoregressive (INAR(1)) model with
zero-inflated (ZI-INAR(1)) and hurdle (H-INAR(1)) innovations is widely
used in studying integer-valued time-series data, such as crime count
and heatwave frequency. This work implemented the INAR(1) models in
Stan.

## Installation

You can install ZIHINAR1 from GitHub with:

``` r
remotes::install_github("fushengyy/ZIHINAR1")
```

Or you can install the released version of ZIHINAR1 from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ZIHINAR1")
```

## Basic Features get_stanfit()

The package contains main function named get_stanfit().

``` r
stan_fit <- get_stanfit(mod_type, distri, y, n_pred = 4,
                        thin = 2, chains = 1, iter = 2000, warmup = iter/2,
                        seed = NA)
```

- mod_type: Character string indicating the model type. Use “zi” for
  zero-inflated models and “h” for hurdle models.

- distri: Character string specifying the distribution. Options are
  “poi” for Poisson or “nb” for Negative Binomial.

- y: A numeric vector of integers representing the observed data.

- n_pred: Integer specifying the number of time points for future
  predictions (default is 4).

- thin: Integer indicating the thinning interval for Stan sampling
  (default is 2).

- chains: Integer specifying the number of Markov chains to run (default
  is 1).

- iter: Integer specifying the total number of iterations per chain
  (default is 2000).

- warmup: Integer specifying the number of warmup iterations per chain
  (default is iter/2).

- seed: Numeric seed for reproducibility (default is NA).

## Example

The following are examples showing how to fit the INAR(1) model when
data is generated from a zero-inflated Negative Binomial (ZINB)
distribution.

``` r
library(ZIHINAR1)
y_data <- data_simu(n = 100, alpha = 0.5, rho = 0.3, theta = c(5, 2), 
                    mod_type = "zi", distri = "nb")
stan_fit <- get_stanfit(mod_type = "zi", distri = "nb", y = y_data, n_pred = 5, 
                        iter = 2000, chains = 1, warmup = 500, 
                        thin = 2, seed = 42)
#> Running /Library/Frameworks/R.framework/Resources/bin/R CMD SHLIB foo.c
#> using C compiler: ‘Apple clang version 17.0.0 (clang-1700.4.4.1)’
#> using SDK: ‘MacOSX26.1.sdk’
#> clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include/"  -I"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppEigen/include/"  -I"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppEigen/include/unsupported"  -I"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/BH/include" -I"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/StanHeaders/include/src/"  -I"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/StanHeaders/include/"  -I"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppParallel/include/"  -I"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c foo.c -o foo.o
#> In file included from <built-in>:1:
#> In file included from /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:22:
#> In file included from /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppEigen/include/Eigen/Dense:1:
#> In file included from /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppEigen/include/Eigen/Core:19:
#> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:679:10: fatal error: 'cmath' file not found
#>   679 | #include <cmath>
#>       |          ^~~~~~~
#> 1 error generated.
#> make: *** [foo.o] Error 1
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000186 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.86 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  501 / 2000 [ 25%]  (Sampling)
#> Chain 1: Iteration:  700 / 2000 [ 35%]  (Sampling)
#> Chain 1: Iteration:  900 / 2000 [ 45%]  (Sampling)
#> Chain 1: Iteration: 1100 / 2000 [ 55%]  (Sampling)
#> Chain 1: Iteration: 1300 / 2000 [ 65%]  (Sampling)
#> Chain 1: Iteration: 1500 / 2000 [ 75%]  (Sampling)
#> Chain 1: Iteration: 1700 / 2000 [ 85%]  (Sampling)
#> Chain 1: Iteration: 1900 / 2000 [ 95%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.663 seconds (Warm-up)
#> Chain 1:                1.76 seconds (Sampling)
#> Chain 1:                2.423 seconds (Total)
#> Chain 1:
get_est(distri = "nb", stan_fit = stan_fit)
#>             Mean         SD    Median       Q2.5      Q97.5      Rhat
#> alpha  0.4411611 0.05248824 0.4436986 0.32853968  0.5305426 0.9987282
#> rho    0.2571094 0.09839077 0.2635160 0.06287134  0.4379304 1.0052384
#> lambda 4.4170520 0.54042828 4.4392172 3.31625159  5.4534335 1.0029409
#> phi    4.5786304 3.44632446 3.6615262 1.33904497 13.6384267 1.0013069
#>        95%_HPD_Lower 95%_HPD_Upper
#> alpha     0.33937746     0.5390050
#> rho       0.05072972     0.4170477
#> lambda    3.33134789     5.4550152
#> phi       0.97038375    10.9572907
get_mod_sel(y = y_data, mod_type = "zi", distri = "nb", stan_fit = stan_fit)
#>       EAIC     EBIC      DIC    WAIC1    WAIC2
#> 1 514.5794 524.9599 518.5852 510.0774 510.4581
get_pred(stan_fit = stan_fit)
#> $summary
#>          Mode Median IQR Min Max
#> y_pred.1    4      5   5   0  19
#> y_pred.2    5      5   5   0  24
#> y_pred.3    6      6   5   0  22
#> y_pred.4    4      5   5   0  24
#> y_pred.5    4      5   5   0  23
#> 
#> $plots
#> $plots[[1]]
```

<img src="man/figures/README-example1-1.png" width="70%" style="display: block; margin: auto;" />

    #> 
    #> $plots[[2]]

<img src="man/figures/README-example1-2.png" width="70%" style="display: block; margin: auto;" />

    #> 
    #> $plots[[3]]

<img src="man/figures/README-example1-3.png" width="70%" style="display: block; margin: auto;" />

    #> 
    #> $plots[[4]]

<img src="man/figures/README-example1-4.png" width="70%" style="display: block; margin: auto;" />

    #> 
    #> $plots[[5]]

<img src="man/figures/README-example1-5.png" width="70%" style="display: block; margin: auto;" />
