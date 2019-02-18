#' @title Koble-Carrigan Isotherm Plot
#' @description Plot of the analysis of Koble-Carrigan Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the plot for koble carrigan isotherm analysis
#' @examples kcarriganplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
kcarriganplot <- function(Ce, Qe){
    x <- 1/ Ce
    y <- 1/Qe
    plot(y~x, xlab="1/ce", ylab="1/qe", main="Koble-Carrigan analysis")
    fit61 <- lm(y~x)
     abline(fit61, col="black")
}
