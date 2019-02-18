#' @title Baudu Isotherm Plot
#' @description Plot of the analysis of Baudu Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the LSRL plot of Baudu isotherm analysis
#'  @export
bauduplot<- function(Ce, Qe){
  x <- 1/Ce
  y <- Qe
  fit10<- lm(y~x)
  plot(y~x, xlab="1/Ce", ylab="Qe", main="Baudu Analysis")
  abline(fit10, col="black")
}
