#' @title Marckzewski Jaroniec Isotherm Plot
#' @description Plot of the analysis of Marckzewski Jaroniec Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the plot for the Marckzewski Jaroniec isotherm analysis
#' @examples marcJplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export

marcJplot <- function(Ce,Qe){
	x <- (1/Ce^2)
y <- Qe

plot(y~x)
fit69 <- lm(y~x)
abline(fit69, main="M-J plot", xlab="1/Ce", ylab="Qe")
}
