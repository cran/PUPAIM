#' @title Sips Isotherm Plot
#' @description Plot of the analysis of Sips Isotherm
#' @param Ce the numerical value for the equilbrium concentration
#' @param Qe the numerical value for the adsorbed concentration
#'
#' @return the LSRL of Sips analysis
#' @examples sipsplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export

sipsplot<- function(Qe,Ce){
	x <- -log10(1/Qe)
	y <- log10(Ce)
	fit77 <- lm(y~x)
plot(y~x, xlab="ln (1/Qe)", ylab="ln Ce", main="Sips Isotherm Plot")
abline(fit77, col="black")
}
