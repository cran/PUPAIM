#' @title Henry Isotherm Plot
#' @description Plot of the analysis of Henry Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the LSRL plot for Henry isotherm analysis
#' @examples henryanalysisplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
henryanalysisplot <- function(Ce,Qe){
	x1 <- Ce
	y1 <- Qe
	plot(x1, y1, main="Henry Isotherm Analysis" , xlab= "Ce" , ylab="Qe", type="p" , col="blue", pch=10)
	fit40 <- lm(y1~x1)
	abline(fit40, col="green")
}
