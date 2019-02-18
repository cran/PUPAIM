#' @title Hill Deboer Isotherm
#' @description Plot of the analysis of Hill Deboer Isotherm
#' @param theta a numeric vector consists of fractional coverage
#' @param Ce a numeric vector consists of equilibrium concentration
#'
#' @return the LSRL plot for Hill Deboer isotherm
#' @examples hdplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
hdplot<- function(theta, Ce){
    x <- theta/(8.314*(273.15+25))
    y <- -log10(Ce)
    fit45 <- lm(y~x)
    plot(y~x, xlab="Theta", ylab="Ce", main="Hill-Deboer Analysis")
abline(fit45 , col="black")
}
