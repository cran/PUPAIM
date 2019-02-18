#' @title BET Isotherm Analysis
#' @description An isotherm that takes account of the possibility that the monolayer in the Langmuir adsorption isotherm can act as a substrate for further adsorption.
#' @param Ce the numerical value for the equilbrium concentration
#' @param Qe the numerical value for the adsorbed concentration
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#' @return the Linear model for the BET isotherm analysis
#'
#' @export
#'
BETanalysis <- function(Ce,Qe){
	x <- Ce
	y <- 1/Qe

	fit11 <- lm(y~x)
	print("BET analysis")
	print(summary(fit11))
	rhs <- function(x,b0,b1){
	((b0-1)/b1*b0)*x
}
	fit12 <- nlsLM(y~ rhs(x,Cs,Qs), start=list(Cs=1, Qs=1))
	print("Parameters")
	print(summary(fit12))


	a <- (summary(fit11))
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
	v <- errors(Ce,Qp)
	print(v)

	}
