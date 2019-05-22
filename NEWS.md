# version 3.0.99 (dev)
- fix bug happening when providing a single Y variable to sglrTheme (issue #1)
- allow to provide folds for cross validation defined by user (issue #2)

# version 3.0
This major version introduces a new feature allowing to group covariates in so called **themes**.

- added `scglrTheme` and `scglrThemeBackward` to handle theme oriented selection
- reworked `multivariateFormula` to handle themes
- added new plots targeting themes
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
