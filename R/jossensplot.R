#' @title Jossens Isotherm Plot
#' @description Plot of the analysis of Jossens Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the LSRL plot for Jossen's isotherm
#' @examples jossensplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
jossensplot<- function(Ce, Qe){
    x <- (Qe)
    y <- log10(Ce/Qe)
    fit49 <- lm(y~x)

    plot(y~x, xlab="Qe", ylab="ln(Ce/Qe)", main="Jossens Analysis")

  abline(fit49, col="blue")



   	}


