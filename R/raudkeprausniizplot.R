#' @title Raudke Prausniiz Isotherm Plot
#' @description Plot of the analysis of Raudke Prausniiz Isotherm
#' @param Ce Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#'
#' @return Provides the plot for Raudke Prausniiz isotherm model
#' @examples radkePplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
radkePplot<- function(Ce,Qe){
	x <- 1/(1+Ce)^2
	y <- 1/Qe
	fit72<- lm(y~x)
plot(y~x, main="Raudke-Prausniiz Analysis", xlab="1/Ce", ylab="1/Qe")
abline(fit72,col="black")
}
