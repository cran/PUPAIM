#' @title Toth Isotherm Plot
#' @description Plot of the analysis of Toth Isotherm
#' @param Ce the numeric value for the equilibrium concentration
#' @param theta the numeric value for the fractional coverage
#'
#' @return the LSRL plot for a toth isotherm analysis
#' @examples tothplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
tothplot <- function(Ce,theta){
    x <- log10(Ce)
    y <- theta
    plot(x, y, xlab="ce", ylab="theta", main="Toth Isotherm Plot")
    fit80 <- lm(y~x)
     abline(fit80, col="black")

}
