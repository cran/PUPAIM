#' @title Fritz Schlunder Isotherm Plot
#' @description Plot of the analysis of Fritz Schlunder Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the LSRL plot fot Fritz Schlunder isotherm analysis
#' @examples fsplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
fsplot <- function(Ce,Qe){

  	x <- Ce

	y <- Qe

	plot(x, y, main ="Fritz-Schlunder Isotherm Plot", xlab= "Ce", ylab= "Qe", type="p", col="blue", pch=5)

	fit31 <- lm(y~x)

	abline(fit31, col="black")


}
