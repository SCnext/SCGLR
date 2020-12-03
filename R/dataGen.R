#' @title Sample dataset of abundance of genera in tropical moist forest
#' @name dataGen
#' @docType data
#' @author CoForChange project
# @keywords datasets
#' @note The use of this dataset for publication must make reference to the CoForChange project.
#' @description dataGen gives the abundance of 8 common tree genera in the tropical moist forest
#' of the Congo-Basin and 58 geo-referenced variables on 2615 8-by-8 km plots
#'  (observations). Each plot's data was obtained by aggregating the data measured on a variable
#'  number of previously sampled 0.5 ha sub-plots. Geo-referenced environmental variables were
#'  used to describe the physical factors as well as vegetation characteristics.
#'  On each plot, 34 physical factors were used pertaining the description of topography, geology, rainfall... 
#'  Vegetation is characterized through 16-days enhanced vegetation index (EVI) data.
#' @references S. Gourlet-Fleury et al. (2009--2014) CoForChange project: \url{http://www.coforchange.eu/}
#' @references C. Garcia et al. (2013--2015) CoForTips project: \url{http://www.cofortips.org/}
#' @format
#' \tabular{ll}{
#'    \code{Y} \tab matrix giving the abundance of 8 common genera (matrix size = 2615*8).\cr
#'    \code{X} \tab matrix of 56 geo-referenced environmental variables (matrix size = 2615*56).\cr
#'    \code{AX} \tab matrix of 2 additionnal explanatory variables (geology and anthropic interference).\cr
#'    \code{offset} \tab sampled area.\cr
#'    \code{random} \tab forest concession id number.\cr
#' }
#'
NULL
