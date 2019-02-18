#' @title Elovich Isotherm Analysis
#' @description The equation that defines this model is based on a kinetic principle which assumes that adsorption sites increase exponentially with adsorption; this implies a multilayer adsorption.
#' @param Qe The numeric value for the equilibrium concentration
#' @param Ce the numeric value for the adsorbed concentration
#'@importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#' @return the regression analysis for the Elovich Isotherm analysis
#' @examples elovichanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
elovichanalysis <- function(Qe,Ce){

	x <- Qe

	y <- log10(Qe/Ce)

	fit16 <- lm(y~x)

	print("Elovich Analysis")

  	print(summary(fit16))

  	a <- (summary(fit16))
  	b <- a$coefficients[2]
  	c <- a$coefficients[1]
  	print(b)
  	print(c)
  	Qp <- function(Ce){
  	  x <- Ce

  	  j<- (b*x) + c

  	  print(j)
  	}
  	b4 <- Qp(Ce)

  	errors <- function(Ce,Qp){

  	  rmse<- rmse(Ce,b4)
  	  mae<- mae(Ce,b4)
  	  mse <- mse(Ce,b4)
  	  rae <-rae(Ce,b4)
  	  colnames(Ce) <- rownames(Ce) <-colnames(Ce)
  	  list("Predicted values", "relative mean squared error" = rmse, "mean absolute error"=mae , "mean squared error"=mse ,  "relative absolute error"=rae)}
  	errors(Ce,Qp)



}
