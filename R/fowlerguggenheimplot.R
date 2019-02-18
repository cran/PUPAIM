#' @title Fowler Guggenheim Isotherm Plot
#' @description Plot of the analysis of Fowler Guggenheim Isotherm
#' @param theta the numerical value for the fractional coverage
#' @param Ce the numerical value for the equilibrium concentration
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the LSRL plot for Flory Huggins isotherm analysis
#' @examples fgplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
fgplot<- function(theta, Ce){
    x <- theta/24.46529
    y <- Ce
    fit25 <- lm(y~x)
    plot(y~x, xlab="Theta", ylab="Ce", main="Fowler-Guggenheim Analysis")
    abline(fit25,col="black")
}
