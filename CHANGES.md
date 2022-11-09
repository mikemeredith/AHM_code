## AHM code

# CHANGES

## 2022-11-09

* In AHM2 chapter 9 sections 1 to 4 and chapter 10 section 2, new scripts are provided which do not require the `RandomFields` package. On Windows, the old `RandomFields` package (version 3.3.14) still works and can be downloaded from [here](https://cran.r-project.org/src/contrib/Archive/RandomFields/). Do not try to use this on a Linux system.


Tested: Windows 10, R 4.2.2, `jagsUI` (1.5.2.9002) and `unmarked` (1.2.5.9004), and up-to-date CRAN versions of other packages.

# 2022-05-04

* Removed all calls to `library(rgdal)` as none were necessary. The `rgdal` package will be retired before the end of 2023; see [here](https://r-spatial.org//r/2022/04/12/evolution.html).

Tested: Windows 10, R 4.3-devel, `jagsUI` (1.5.2), and up-to-date CRAN versions of other packages.

## 2021-07-12

* Add script to plot the figures in AHM2 section 8.2.3 (Figs 8.4 to 8.7) with credible intervals.

## 2021-06-25

* Incorporate [errata](https://sites.google.com/site/appliedhierarchicalmodeling/errata) into code for sections AHM2 3.4.4 and 3.4.5.

## 2021-06-21

* Changes to calls to `jagsUI::traceplot` to accommodate the new version (1.5.2).

`traceplot` and `densityplot` now plot multiple nodes in the same window, so calling `par(mfrow=...)` before calling the plotting function is no longer necessary or effective. The default with > 4 nodes to plot is a 3 x 3 layout. Other layouts can be specified with the `layout` argument.

Tested: Windows 10, R 4.2-devel, and up-to-date CRAN versions of other packages, including`jagsUI` (1.5.2).

## 2021-05-23

* Minor corrections to several files, mostly typos.
* The testing routine now compares the values produced with the output from a prior run. Previously it only checked for errors in execution.

Tested: Windows 10, R 4.2-devel, GitHub version `jagsUI` (1.5.1.9102) and up-to-date CRAN versions of other packages, including `unmarked` (1.1.1) and `AHMbook` (0.2.3).

## 2021-05-14

* Fixed incorrect description of `lambda` in comments in AHM2_02.05.1+2.R (thanks to José Jiménez).
* Fixed run time for the same file.

## 2021-04-28

* Fixed restoration of plotting parameters after calls to `par` in multiple files, as autocheck now reports when this is not done.

Tested: Windows 10, R 4.1.0-alpha, GitHub versions of `AHMbook` (0.2.2.9001), `jagsUI` (1.5.1.9101) and `unmarked` (1.0.1.9014), up-to-date CRAN versions of other packages.

## 2021-02-24 to 27

* Added scripts with NIMBLE code for CAR models in AHM2 sections 3.4.4, 3.4.5, 9.4.1, 9.4.3 and 09.5.

## 2021-02-19

* Removed all non-ASCII characters (mostly smart quotes) from _comments_ in BUGS/JAGS code. These did not affect running the models, but prevented reloading the output saved in `.RData` files.

## 2021-02-18

* Fixed number of cores used by `AICcmodavg::Nmix.gof.test` and `unmarked::parboot` when `parallel=TRUE`; the default is to use all-but-one of the cores on the machine, and can crash when other applications are running.

Tested: Windows 10, R 4.0.4, GitHub versions of `AHMbook` (0.2.2.9001), `jagsUI` (1.5.1.9101) and `unmarked` (1.0.1.9011), up-to-date CRAN versions of other packages.

## 2021-02-13

* Added bonus script to run the simulation in AHM1 10.7 in JAGS

## 2020-10-08 to 27

* Added code for AHM2 Figures in chapters 2, 7, 11
* Added simulation code for AHM2 1.7.1
* Removed unnecessary call to `library(coda)` in AHM1 11.7 to avoid clashes with `jagsUI::traceplot`.
* Tidied up plotting code in AHM1 11.7.

Tested: Windows 10, R 4.0.3, jagUI 1.5.1 patched to `ask` only for interactive devices, up-to-date CRAN versions of other packages.

## 2020-07-14 to 31

* Added AHM2 chapters 7 to 11.

## 2020-06-27

Many small changes to AHM1 scripts

* Added code execution times to scripts which take more than a minute or so to run.
* Added `par(op)` to reset plotting parameters after each block of plotting code.
* Changed indenting in many places to make it consistent across scripts.
* Broke up long lines to keep line length <= 80 characters as far as possible.
* `T` and `F` replaced with `TRUE` and `FALSE` wherever appropriate.

## 2020-06-25

* Added TO_DO file.
* Added AHM2 chapters 5 and 6, including bonus code and figure code.

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
* Plotting calls with `ask` argument get `dev.interactive()`, this means R only asks for page confirmation for screen devices. (Testing sends plots to a .PDF file.)
* In a few cases, calls to `jagsUI::jags` have been changed to `parallel = TRUE` to save time.
* Numbers of iterations have been reduced for some simulations and MCMC runs that were taking several hours.

Tested: Windows 10, R 3.6.2, jagUI 1.5.1 patched to `ask` only for interactive devices, up-to-date CRAN versions of other packages.

# 2019-08-09

Changes to scripts so that each script can be run in a new R session. This involved inserting additional code from earlier sections at the top of the script or loading saved output from previous sections.

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

