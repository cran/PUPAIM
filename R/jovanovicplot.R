#' @title Jovanovic Isotherm Plot
#' @description Plot of the analysis of Jovanovic Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the LSRL plot for Jovanovic isotherm analysis
#' @examples jovanovicplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
jovanovicplot <- function(Ce,Qe)
  {
	x5<-Ce
	y5<-log10(Qe)
	plot(x5, y5, main = "Jovanovic Isotherm", xlab = "Ce", ylab= "ln(Qe)", type = "p", col ="blue", pch =10)
	fit52 <- lm(y5~x5)
	abline(fit52, col = "green")
  }

