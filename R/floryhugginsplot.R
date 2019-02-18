#' @title Flory Huggins Isotherm Plot
#' @description Plot of the analysis of Flory Huggins
#' @param Ce the numerical value for the equilibrium capacity
#' @param theta the numerical value for the fractional coverage
#'
#' @return the LSRL plot for Flory Huggins isotherm analysis
#'
#' @export
#'
fhplot<- function(Ce, theta){
    x <- log10(1-theta)
    y <- log10(theta/Ce)
    fit22 <- lm(y~x)
    plot(y~x, xlab="ce", ylab="theta", main="Flory Huggins Analysis")
    abline(fit22, col="black")
}
