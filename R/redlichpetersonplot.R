#' @title Redlich Peterson Isotherm Plot
#' @description Plot of the analysis of Redlich Peterson Isotherm
#' @param Ce Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#'
#' @return Provides the plot for Redlich Peterson isotherm model
#' @examples rpplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
rpplot<- function(Ce, Qe){
    x <- log10(Ce)
    y <- Qe
    fit75 <- lm(y~x)
    plot(y~x, xlab="Ce", ylab="Qe", main="Redlich-Peterson Analysis")
    abline(fit75, col="black")

}
