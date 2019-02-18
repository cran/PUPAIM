#' @title Elovich Isotherm Plot
#' @description Plot of the analysis of Elovich Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the LSRL plot for Elovich isotherm analysis
#' @examples elovichplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
elovichplot <- function(Qe,Ce){

  	x <- Qe

	y <- log10(Qe/Ce)

	plot(x, y, main ="Elovich Isotherm Plot", xlab= "Ce", ylab= "log10(Qe/Ce)", type="p", col="blue", pch=10)

	fit19 <- lm(y~x)

	abline(fit19, col="black")


}
