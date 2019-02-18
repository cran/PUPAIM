#' @title BET Isotherm Plot
#' @description Plot of the analysis of BET Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#
#' @return the plot for the BET adsorption isotherm
#' @examples betplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'

betplot <- function(Ce,Qe){
  x <- Ce
  y <- 1/Qe

  fit84 <- lm(y~x)
  plot(y~x, main="BET Isotherm Analysis", xlab= Ce, ylab = "1/Qe")
  abline(fit84)
}
