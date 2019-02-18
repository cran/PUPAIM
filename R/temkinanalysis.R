#' @title Temkin Isotherm Analysis
#' @description takes into account the effects of indirect adsorbate/adsorbate interaction on the adsorption process
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the linear regression and the parameters for the Temkin Isotherm Analysis
#' @examples temkinanalysis (c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
temkinanalysis <- function(Ce,Qe){
   x <- log(Ce)
   y <- Qe

   fit <- lm(y~x)
   print("Temkin Analysis")
   print(summary(fit))

   rhs <- function(x, b0,b1){
     b1*(log(b0)) + b1*x
   }
   fit1 <- nlsLM(y~rhs(x, Kt,b), start=list(Kt=1,b=1), trace = TRUE)
   print("Parameters")
   print(summary(fit1))
}
