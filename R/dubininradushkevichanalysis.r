#' @title Dubinin Radushkevich Isotherm Analysis
#' @description It is only suitable for intermediate range of adsorbate concentrations because it exhibits unrealistic asymptotic behavior and does not predict Henry’s laws at low pressure.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the linear regression and the parameters for the Dubinin Radushkevich isotherm analysis
#' @examples dubininradanalysis (c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
dubininradanalysis <- function(Ce,Qe){
	x <- Ce^2
	y <- log10(Qe)
	fit13 <- lm(y~x)
	print("Dubinin-Radushkevich Analysis")
	print(summary(fit13))

	rhs <- function(x, b0,b1){
	log10(b1) - b0*x
	}
	fit14 <- nlsLM(y~ rhs(x, beta, qm), start=list(beta=3, qm=4), trace = TRUE)
	print("Parameters")
	print(summary(fit14))

	a <- (summary(fit13))
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
