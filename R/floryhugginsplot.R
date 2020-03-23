#' @title Flory Huggins Isotherm Plot-deprecated
#' @description Plot of the analysis of Flory Huggins
#' @param Ce the numerical value for the equilibrium capacity
#' @param theta the numerical value for the fractional coverage
#' @return the LSRL plot for Flory Huggins isotherm analysis
#'
#' @name fhplot-deprecated
#' @usage fhplot(Ce, theta)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{fhplot}:
#' For \code{fhplot}, use \code{\link{fhanalysis}}.
#'
#' @export
fhplot<- function(Ce, theta){
  .Deprecated("fhanalysis")
  x <- log10(1-theta)
  y <- log10(theta/Ce)
  fit22 <- lm(y~x)
  plot(y~x, xlab="ce", ylab="theta", main="Flory Huggins Analysis")
  abline(fit22, col="black")
}
