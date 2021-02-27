# AHM code

The two-volume work *Applied Hierarchical Modeling in Ecology: Analysis of distribution, abundance and species richness in R and BUGS* by Marc Kéry and Andy Royle contains lots of R and BUGS code.

The R package `AHMbook`, available on CRAN, has all the data sets and the custom functions used in the books. Commented code for the functions is on GitHub [here](https://github.com/mikemeredith/AHMbook).

This repository has all the code in the printed books, plus code referred to as "available on the website" but not printed. The aim is to have code which works with current versions of R, JAGS and contributed R packages. The code is regularly tested and updated code inserted, with the original printed code retained but commented out with `#`. Please open an issue if you find other code which does not work.

In addition to these updates, some code has been inserted:
* Code added at the top of the script to recreate or reload objects from previous sections; each script is self-contained.
* Some changes have been made to facilitate automated checking of scripts, in particular reductions in the number of iterations for simulations, bootstraps and MCMC runs.
* After long runs of `unmarked`, `JAGS` or `WinBUGS`, I've inserted code to save the results to `RData` files.

Additional code and comments are marked off with twiddly lines like this:
```
#~~~~ oldfunction has been replaced with newfunction ~~~~~~~
# oldfunction(foo)
newfunction(foo)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

## Avoiding WinBUGS

WinBUGS is not essential to work through the code. In most cases, `jagsUI::jags` is a drop-in replacement for `R2WinBUGS::bugs`. JAGS does not do spatial autocorrelation (CAR) models, as used in AHM2 chapters 3 and 9; for those, the `nimble` package can be used, and alternative scripts are provided.

## Volume 1 (AHM1)

The code for the first volume (AHM1), with updates up to 2017, is available as a single huge text file on the [main book web page](http://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/).

## Volume 2 (AHM2)

The book appeared in October 2020 (with copyright dated 2021). The code here is based on the final proofs.
