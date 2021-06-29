# version 3.0.9000 (dev)
- fix bug happening when providing a single Y variable to `sglrTheme` (issue #1)
- allow to provide folds for cross validation defined by user (issue #2)
- individual beta and gamma for themes were not properly computed (issue #5) 
- improve evaluation of linear predictors in `theme` function (issue #8)
- remove dependency on `expm` package (issue #10)
- error message for lack of convergence were improved and `scglrCrossVal` no longer
  stops when a fold does not converge (issue #12)
- add a plot for `scglrCrossVal`
- fix contrasts error in `scglrTheme` (issue #20) 

## New features
- preliminary integration of code from `SCnext/mixedSCGLR` written by Jocelyn Chauvet (issue #11)
- moved parallelization to `future` framework. Need `future.apply` package (issue #16)
- added progress to `scglrCrossVal` and `scglrThemeBackward` functions using `progressr` package (issue #17)

## Deprecations
- parameter `nfolds` of `scglrCrossVal` as been renamed to `folds` and now accepts also vectors.
For now a warning is issued and provided value is used but update your code as it will be removed
in a next version.
- parameter `mc.cores` of `scglrCrossVal` is now deprecated as of the move to `future` framework for parallelization.

# version 3.0
This major version introduces a new feature allowing to group covariates in so called **themes**.

- added `scglrTheme` and `scglrThemeBackward` to handle theme oriented selection
- reworked `multivariateFormula` to handle themes
- added new plots targeting themes

## Deprecations
- deprecated `barplot` in favor of `screeplot` (same parameters)

# version 2.1
- Removed LPLS legacy method
- changed from ING to PING

# version 2.0
- New method is available : SR (Structural Relevance) see vignette
- Major rewrite of plot styling (not backward compatible)
- Various fixes and improvements (especially when dealing with
a single dependant variable)

# version 1.0
Initial version of SCGLR

- method LPLS (Local PLS)
