% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.r
\name{plot.SCGLR}
\alias{plot.SCGLR}
\title{SCGLR generic plot}
\usage{
\method{plot}{SCGLR}(x, ..., style = getOption("plot.SCGLR"), plane = c(1, 2))
}
\arguments{
\item{x}{an object from SCGLR class.}

\item{\dots}{optional arguments (see \link{customize}).}

\item{style}{named list of values used to customize the plot (see \link{customize})}

\item{plane}{a size-2 vector (or string with separator) indicating which components 
are plotted (eg: c(1,2) or "1,2" or "1/2").}
}
\value{
an object of class \code{ggplot}.
}
\description{
SCGLR generic plot
}
\examples{
\dontrun{
library(SCGLR)

# load sample data
data(genus)

# get variable names from dataset
n <- names(genus)
ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
nx <- n[-grep("^gen",n)]   # X <- remaining names

# remove "geology" and "surface" from nx
# as surface is offset and we want to use geology as additional covariate
nx <-nx[!nx\%in\%c("geology","surface")]

# build multivariate formula
# we also add "lat*lon" as computed covariate
form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))

# define family
fam <- rep("poisson",length(ny))

genus.scglr <- scglr(formula=form,data = genus,family=fam, K=4,
 offset=genus$surface)

summary(genus.scglr)

barplot(genus.scglr)

plot(genus.scglr)

plot(genus.scglr, predictors=TRUE, factor=TRUE)

pairs(genus.scglr)

}
}
