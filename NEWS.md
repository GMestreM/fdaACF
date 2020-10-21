# fdaACF 1.0.1

Fix minor bugs:

* Fix the automatic selection of variable `n_harm` in function `FTS_identification`. The previous version obtained the first 10 principal components of the functional time series and looked for the number of components neccesary to explain more than 95 % of the variance of the data. This approach presents several drawbacks: perhaps more than 10 PC are neccesary to explain more than 95% of the variance, or the series could have less than 10 discretization points.
The updated version of function `FTS_identification` now uses the length of the input variable `v` to decide the number of PC to try and gives a warning message if 95 % of the variance cannot be explained by the maximum number of PC tried.

# fdaACF 1.0.0

New functionalities added in this release:

* Add function FTS_identification. This function is a wrapper for functions obtain_FACF and obtain_FPACF, and shows both the FACF and the FPACF of a given functional time series.
* Add citation information for the paper "Functional time series model identification and diagnosis by means of auto- and partial autocorrelation analysis" [<doi:10.1016/j.csda.2020.107108>](https://doi.org/10.1016/j.csda.2020.107108). This paper contains the theoretical foundations of the autocorrelation functions implemented in this package.
* Fix minor bugs present in version 0.2.0.

# fdaACF 0.2.0

New functionalities added in this release:

* Add function obtain_FPACF. This allows the user to estimate the functional partial autocorrelation function of a given functional time series.


# fdaACF 0.1.0

Initial release of the fdaACF library to CRAN.

This release contains the core functionalities of the fdaACF library.
Available at: https://CRAN.R-project.org/package=fdaACF




