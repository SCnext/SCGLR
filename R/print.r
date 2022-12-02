#' @export
#' @method print SCGLR
#' @keywords internal
#' @title Print SCGLR object
#' @description Prints inertia per component and residual and null deviance for each Y.
#' @param x object of class 'SCGLR', usually a result of running \code{\link{scglr}}.
#' @param \dots Not used.
print.SCGLR <- function(x, ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "","\n")
  cat("\nInertia:\n")
  print.default(x$inertia,print.gap=2)
  cat(sprintf("\nNull deviance on %d degrees of freedom:\n",attr(x$deviance.null,"df")))
  attributes(x$deviance.null) <- NULL
  print.default(x$deviance.null,print.gap=2)
  cat(sprintf("\nResidual deviance on %d degrees of freedom:\n",attr(x$deviance.residual,"df")))
  attributes(x$deviance.residual) <- NULL
  print.default(x$deviance.residual,print.gap=2)
  invisible(x)
}
