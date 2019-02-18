#' @title Temkin Isotherm Analysis
#' @description takes into account the effects of indirect adsorbate/adsorbate interaction on the adsorption process
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the plot for the Temkin Isotherm Analysis
#' @examples temkinplot (c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
temkinplot <- function(Ce,Qe){
  x <- log(Ce)
  y <- Qe
  plot(y~x,main="Temkin plot", xlab="Ce", ylab="Qe")
  fit <- lm(y~x)
  abline(fit, col="black")
}
