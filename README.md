<!-- File generated from README.Rmd. Changes must be done from there -->

# SCGLR

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/SCGLR)](https://cran.r-project.org/package=SCGLR)

## Introduction

**SCGLR** is an open source implementation of the Supervised Component
Generalized Linear Regression (Bry et al. [2013](#ref-bry13),
[2016](#ref-bry16), [2018](#ref-bry18)), which identifies, among a large
set of potentially multicolinear predictors, the strong dimensions most
predictive of a set of responses.

**SCGLR** is an extension of partial least square regression (PLSR) to
the uni- and multivariate generalized linear framework. PLSR is
particularly well suited for analyzing a large array of explanatory
variables and many studies have demonstrated its predictive performance
in various biological fields such as genetics (Boulesteix and Strimmer
[2007](#ref-boulesteix07)) or ecology (Carrascal, Galván, and Gordo
[2009](#ref-carrascal09)). While PLSR is well adapted for continuous
variables, maximizing the covariance between linear combination of
dependent variables, and linear combinations of covariates, **SCGLR** is
suited for non-Gaussian outcomes and non-continuous covariates.

**SCGLR** is a model-based approach that extends PLS (Tenenhaus
[1998](#ref-tenenhaus98)), PCA on instrumental variables (Sabatier,
Lebreton, and Chessel [1989](#ref-sabatier89)), canonical correspondence
analysis (Ter Braak [1987](#ref-terbraak87)), and other related
empirical methods, by capturing the trade-off between goodness-of-fit
and common structural relevance of explanatory components. The notion of
structural relevance has been introduced (Bry and Verron
[2015](#ref-bry15)).

**SCGLR** can deal with covariates partitioned in several groups called
“themes”, plus a group of additional covariates. Each theme is
searched for orthogonal components representing its variables in the
model, whereas the additional covariates appear directly in the model,
without the mediation of a component (Bry et al. [2019](#ref-bry19)).

## Installation

``` r
# Install release version from CRAN
install.packages("SCGLR")

# Install development version from GitHub
devtools::install_github("SCnext/SCGLR")
```

## Main functions and works in progress

**SCGLR** is designed to deal with outcomes from multiple distributions:
Gaussian, Bernoulli, binomial and Poisson separately or simultaneously
(Bry et al. [2013](#ref-bry13)). Moreover **SCGLR** is also able to deal
with multiple conceptually homogeneous explanatory variable groups (Bry
et al. [2018](#ref-bry18)).

**SCGLR** is a set of **R** functions illustrated on a floristic data
set, *genus*. `scglr` and `scglrTheme` are respectively dedicated to
fitting the model with one or more thematic group of regressors.
`scglrCrossVal` and `scglrThemeBackward` are respectively dedicated to
selecting the number of components. `print`, `summary` and `plot`
methods are also available for the `scglr` and `scglrTheme` function
results.

Different works are in progress both dealing for instance with the
inclusion of random effects extending **SCGLR** to the generalized
linear mixed model framework (Chauvet, Trottier, and Bry
[2018](#ref-chauvet18)[a](#ref-chauvet18),
[2018](#ref-chauvet18b)[b](#ref-chauvet18b)), or the Cox regression
model.

## References

<div id="refs" class="references">

<div id="ref-boulesteix07">

Boulesteix, Anne-Laure, and Korbinian Strimmer. 2007. “Partial Least
Squares: A Versatile Tool for the Analysis of High-Dimensional Genomic
Data.” *Briefings in Bioinformatics* 8 (1): 32–44.
<http://bib.oxfordjournals.org/content/8/1/32.short>.

</div>

<div id="ref-bry19">

Bry, Xavier, Catherine Trottier, Frédéric Mortier, and Guillaume Cornu.
2019. “Component-Based Regularisation of a Multivariate GLM with a
Thematic Partitioning of the Explanatory Variables.” *Statistical
Modelling* 19 (0): 00–00 (to appear). <https://doi.org/TO BE ADDED>.

</div>

<div id="ref-bry18">

Bry, X., C. Trottier, F. Mortier, and G Cornu. 2018. “Component-Based
Regularisation of a Multivariate Glm with a Thematic Partitioning of the
Explanatory Variables.” *Statistical Modelling*, In press.

</div>

<div id="ref-bry16">

Bry, X., C. Trottier, F. Mortier, G. Cornu, and Verron T. 2016.
“Supervised-Component-Based Generalised Linear Regression with
Multiple Explanatory Blocks: THEME-Scglr.” In *The Multiple Facets of
Partial Least Squares and Related Methods*, edited by H. Abdi, V.E.
Vinzi, V. Russolillo, G. Saporta, and L Trinchera, 141–54. Switzerland:
Springer Proceedings in Mathematics & Statistics.

</div>

<div id="ref-bry13">

Bry, X., C. Trottier, T. Verron, and F. Mortier. 2013. “Supervised
Component Generalized Linear Regression Using a Pls-Extension of the
Fisher Scoring Algorithm.” *Journal of Multivariate Analysis* 119:
47–60.
<http://www.sciencedirect.com/science/article/pii/S0047259X13000407>.

</div>

<div id="ref-bry15">

Bry, X., and T Verron. 2015. “THEME: THEmatic Model Exploration Through
Multiple Co-Structure Maximization.” *Journal of Chemometrics* 29 (12):
637–47. <http://onlinelibrary.wiley.com/doi/10.1002/cem.2759/full>.

</div>

<div id="ref-carrascal09">

Carrascal, Luis M., Ismael Galván, and Oscar Gordo. 2009. “Partial Least
Squares Regression as an Alternative to Current Regression Methods Used
in Ecology.” *Oikos* 118 (5): 681–90.
<http://onlinelibrary.wiley.com/doi/10.1111/j.1600-0706.2008.16881.x/full>.

</div>

<div id="ref-chauvet18">

Chauvet, J., C. Trottier, and X Bry. 2018a. “Component-Based
Regularisation of Multivariate Generalised Linear Mixed Models.”
*Journal of Computational and Graphical Statistics*, In press.

</div>

<div id="ref-chauvet18b">

———. 2018b. “Regularisation of Generalised Linear Mixed Models with
Autoregressive Random Effect.” *Journal of Computational and Graphical
Statistics*, In prep.

</div>

<div id="ref-sabatier89">

Sabatier, R., J. D. Lebreton, and D. Chessel. 1989. “Principal Component
Analysis with Instrumental Variables as a Tool for Modelling Composition
Data.” *Multiway Data Analysis*, 341–52.

</div>

<div id="ref-tenenhaus98">

Tenenhaus, M. 1998. *La Régression PLS: Théorie et Pratique*. Paris:
Editions Technip.
<https://books.google.fr/books?hl=fr&lr=&id=OesjK2KZhsAC&oi=fnd&pg=PA1&dq=Tenenhaus+PLS&ots=EvUst85CEP&sig=EpksVNlZFUVoYLX7JX952PIGaHU>.

</div>

<div id="ref-terbraak87">

Ter Braak, Cajo JF. 1987. “The Analysis of Vegetation-Environment
Relationships by Canonical Correspondence Analysis.” In *Theory and
Models in Vegetation Science*, 69–77. Springer.
<https://link.springer.com/chapter/10.1007/978-94-009-4061-1_7>.

</div>

</div>
