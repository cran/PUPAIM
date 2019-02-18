#' @title Halsey Isotherm Plot
#' @description Plot of the analysis of Halsey Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the LSRL plot for Halsey isotherm analysis
#' @examples halseyplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
halseyplot <- function(Ce,Qe){
    x<-log10(Ce)
    y <- Qe
    plot(x, y, xlab="ln ce", ylab="Qe", main="Halsey Analysis")
    fit34<- lm(y~x)
     abline(fit34, col="black")


}
