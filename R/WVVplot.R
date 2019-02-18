#' @title Weber Van Vliet Isotherm Plot
#' @description Plot of the abalysis of Weber Van Vliet Isotherm
#'
#' @param Qe the numeric value for the adsorbed concentration
#' @param Ce the numeric value for the equilibrium concentration
#'
#' @return the linear plot for the Weber Van Vliet Isotherm
#' @examples WVVplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
WVVplot <- function(Qe, Ce){
 x <- Qe
	y <- Ce

plot(y~x, main="WVV analysis", xlab="Ce", ylab="Qe")
fit83 <- lm(y~x)
abline(fit83, col="black")
}
