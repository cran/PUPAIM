#' @title Fractional Power Kinetic Model
#' @description Calculates Kt parameter of the Fractional power kinetic model.
#' @param t duration of the experiment
#' @param qt the numerical value for the equilibrium capacity
#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @importFrom minpack.lm "nlsLM"
#' @return the regression analysis for the first order kinetics
#' @examples fractionalpower (c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
fractionalpower <- function(t, qt){
  x <- log10(t)
  y <- log10(qt)

  fit <- lm(y~x)
  print(summary(fit))

  rhs <- function(x,b0,b1){
    log10(b0) + b1*x
  }
   fit1 <- nlsLM(y~rhs(x,a,b), start=list(a=1,b=1))
   print("parameters")
   print(summary(fit1))
}
