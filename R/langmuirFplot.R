#' @title Langmuir-Freundlich Isotherm Plot
#' @description Plot of the analysis of Langmuir-FreundLich Isotherm
#' @param Ce the numerical value for the equilbrium concentration
#' @param Qe the numerical value for the adsorbed concentration
#'
#' @return the LSRL for the Langmuir-Freundlich isotherm analysis
#' @examples langmuirFplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
langmuirFplot <- function(Ce,Qe){
	x <- 1/(Ce)^2
	y <- Qe
	fit66 <- lm(y~x)
	plot(y~x, main="Langmuir-Freundlich Analysis", xlab="1/Ce^2", ylab = "Qe")
	abline(fit66,col="black")
}
