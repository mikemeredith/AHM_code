## AHM code

# CHANGES

## 2020-06-21

More material for AHM2

* Added plotting code for most of the figures in chapters 1 to 3.
* Added chapter 4, including bonus code and figure code.

## 2020-06-01

Testing with R 4.0.1 RC and latest versions of packages.

* Change of default for `stringsAsFactors` in `data.frame`: `stringsAsFactors=TRUE` needed in some places.

## 2020-02-06

Preliminary code for AHM2 chapters 1 to 3 added. These all work properly but may need some tidying up.

## 2019-12-14

A slew of changes to enable scripts to run reasonably quickly and without human intervention for testing purposes. The main features are:

* `R2WinBUGS::bugs` calls now all have `debug = FALSE`, original lines with `debug = TRUE` commented out.
* Calls to `browser()` commented out, replaced with a call to `devAskNewPage()`.
* Plotting calls with `ask` argument get `dev.interactive()`, this means R only asks for page confirmation for screen devices.
* In a few cases, calls to `jagsUI::jags` have been changed to `parallel = TRUE` to save time.
* Numbers of iterations have been reduced for some simulations and MCMC runs that were taking several hours.

Tested: Windows 10, R 3.6.2, jagUI 1.5.1 patched to `ask` only for interactive devices, up-to-date CRAN versions of other packages.

# 2019-08-09

Changes to scripts so that each script can be run in a new R session. This involved inserting additional code from earlier section at the top of the script or loading saved output from previous sections.

In checking the scripts, I found some of the original code no longer worked and provided new code that does work. For example:

* `AICcmodavg::Nmix.gof.test` and `unmarked::parboot` now have an argument `parallel` with default `TRUE`, but some original code will not run in parallel.
* The output object for `AICcmodavg::predict` has changed.
* From R 3.6.0 a new random number generator was introduced as the default for `sample` and relatives.
* In sections 11.7-11.9 the order of the columns in `all10` is now different, so column numbers are now wrong.
* Shape files for the plots of Swiss maps referred to in the code have not been made available.

Tested: Windows 10, R 3.6.1 and up-to-date CRAN versions of packages.

# 2019-08-05

Uploaded Vol 1 code from 2017-05-19, chopped into scripts for each main section with Errata up to 2019-08-05.

# 2019-08-02

Test commit with code for Vol 1 chapter 1 and README.md file.

