#' @title Temkin Isotherm Analysis Plot-deprecated
#' @description takes into account the effects of indirect adsorbate/adsorbate interaction on the adsorption process
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom minpack.lm "nlsLM"
#' @return the plot for the Temkin Isotherm Analysis
#'
#' @name temkinplot-deprecated
#' @usage temkinplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL
#' @rdname PUPAIM-deprecated
#' @section \code{temkinplot}:
#' For \code{temkinplot }, use \code{\link{temkinanalysis}}.
#'
#' @export
temkinplot <- function(Ce,Qe){
  .Deprecated("temkinanalysis")
  x <- log(Ce)
  y <- Qe
  plot(y~x,main="Temkin plot", xlab="Ce", ylab="Qe")
  fit <- lm(y~x)
  abline(fit, col="black")
}
