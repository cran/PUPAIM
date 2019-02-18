#' @title Kiselev Isotherm Plot
#' @description Plot of the analysis of Kiselev Isotherm
#' @param Kh the numerical value for the kiselev constant
#' @param Ce the numerical value for the equilibrium capacity
#'
#' @return the plot for kiselev isotherm analysis
#' @examples kiselevplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
kiselevplot <- function(Kh,Ce){
  x<- Kh

  y<-1/Ce
     plot(x, y, xlab="Kh", ylab="1/Ce", main="Kiselev Analysis")
    fit58 <- lm(y~x)
     abline(fit58, col="black")

}
