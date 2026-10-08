
# popClumpTestR

<!-- badges: start -->
<!-- badges: end -->

The goal of popClumpTestR is to ...

## Installation

You can install the development version of popClumpTestR from [GitHub](https://github.com/) with:

``` r
devtools::install_github("Orrslab/popClumpTestR", build_vignettes = TRUE)
library(popClumpTestR)
```

Enables simulating X-Y population point data to test robustness of methods to recognize clumpedness. Given Quantification Methods (QMs), it is possible to simulate many populations at different levels of clumping. Afterwards, it is possible to generate Large Random-like Population Distributions (LRPDs) of your QMs - a baseline for random populations of 100 individuals in high spatial randomness, enabling a distribution per CM of random values. A high percentile per CM is marked as a Minimum Non-Random Threshold (MNRT), and a Non-Random Percent (NRP) of the simulation data is then calculated. Display functions (using ggplot2) then enable to qualitatvely ascertain which CMs are good enough, with differentiability of NRPs at different levels of clumping, and also to ascertain a minimum sample size.
Note that the locations of simulations and example data are in Israel Transverse Mercator (ITM) - https://en.wikipedia.org/wiki/Israeli_Transverse_Mercator. This is to simulate localized coordinates - but technically any localized coordinate system will do.

For in-depth description and explanations, please see the introduction.Rmd.

