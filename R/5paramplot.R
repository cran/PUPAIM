#' @title Five Parameter Isotherm Plot
#' @description Plot of the analysis of Five Parameter Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'@importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#'
#'@return the linear regression and the parameters for the Five parameter isotherm analysis
#' @examples fivePplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
fivePplot <- function(Ce,Qe){
	x <- Ce^2
	y <- Qe
  fit7 <- lm(y~x)
	plot(y~x, main="Five parameter Analysis", xlab="Ce^2",ylab = "Qe")
	abline(fit7, col="black")
}
