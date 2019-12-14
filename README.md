# AHM code

The book *Applied Hierarchical Modeling in Ecology: Analysis of distribution, abundance and species richness in R and BUGS* by Marc Kéry and Andy Royle contains lots of R and BUGS code. The code for the first volume (AHM1), with some updates, is available as a single huge text file on the [main book web page](http://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/). Volume 2 (AHM2) will shortly appear with more code, and the plan is to make this available plus additional code not included in the book.

The initial post was the 2017 code, chopped into separate scripts for each main section of the book. See the CHANGES file for details of changes since then. The main changes are:

* Code added at the top of the script to recreate or reload objects from previous sections.
* Some original code no longer works with current versions of R and packages. New code which does work has been inserted and the old code commented out with `#`.
* Some changes have been made to facilitate automated checking of scripts.
* After long runs of `unmarked`, `JAGS` or `WinBUGS`, I've inserted code to save the results to `RData` files.

Additional comments and new code are marked off with twiddly lines like this:
```
# ~~~~ oldfunction has been replaced with newfunction ~~~~~~~
# oldfunction(foo)
newfunction(foo)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

Please open an issue if you find other code which does not work.
