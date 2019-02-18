#' @title Halsey Isotherm
#' @description Used to evaluate multilayer adsorption at a relatively large distance from the surface.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Halsey isotherm analysis
#'
#' @export
#'
halseyanalysis <-function(Ce,Qe){

	x1 <- log10(Ce)

	y2 <- Qe

	fit32 <- lm(y2~x1)
	print("Halsey Analysis")
	print(summary(fit32))

	rhs <- function(x, b0,b1){

	  ((1/b0)*b1 - (1/b0*x))

}

	fit33 <- nlsLM(y2~rhs(x1,nh,Kh), start=list(nh= 75709.31,Kh=-49715.73), trace=TRUE)
print("Parameters")
	print(summary(fit33))

	a <- (summary(fit32))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit33))
	e <- d$coefficients[1]

	print(b)
	print(c)
	Qp <- function(Qt){
	  x <- Qt

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
	deltaG(25,e)

}


