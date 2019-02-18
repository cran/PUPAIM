#' @title Hill Isotherm Plot
#' @description Plot of the analysis of Hill Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the LSRL plot for Hill isotherm analysis
#' @examples hillplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
hillplot<- function(Ce,Qe){
    x <- log(Ce)
    y <- log(Qe/.5-Qe)
    fit46 <- lm(y~x)
    plot(y~x, xlab="log(Ce)", ylab="Qe", main="Hill Analysis")
    abline(fit46,col="black")

}
