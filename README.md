# AHM code

The book *Applied Hierarchical Modeling in Ecology: Analysis of distribution, abundance and species richness in R and BUGS* by Marc Kéry and Andy Royle contains lots of R and BUGS code. The code for the first volume (AHM1), with some updates, is available as a single huge text file on the [main book web page](http://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/). Volume 2 (AHM2) will shortly appear with more code, and the plan is to make this available plus additional code not included in the book.

The initial post was the 2017 code, chopped into separate scrips for each main section of the book. Scripts have since been updated so that they work: you should be able to open a new R session and `source` a script. That means adding at least calls to `library` at the top of the script and often recapping code from other sections to create data objects.

In checking the scripts (on Windows 10 with R 3.6.1 and up-to-date packages), I found some of the original code no longer worked. For example:

* `AICcmodavg::Nmix.gof.test` and `unmarked::parboot` now have an argument `parallel` with default `TRUE`, but some original code will not run in parallel.
* The output object for `AICcmodavg::predict` has changed.
* From R 3.6.0 a new random number generator was introduced as the default for `sample` and relatives.
* In sections 11.7-11.9 the order of the columns in `all10` is now different, so column numbers are now wrong.
* Shape files for the plots of Swiss maps referred to in the code have not been made available.

After long runs of JAGS or WinBUGS, I've inserted code to save the results to RData files.

Additional comments and new code are marked off with twiddly lines like this:
```
# ~~~~ oldfunction has been replaced with newfunction ~~~~~~~
# oldfunction(foo)
newfunction(foo)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

Please open an issue if you find other code which does not work.
