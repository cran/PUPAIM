#' @title  Weber Van Vliet Isotherm Analysis
#' @description An empirical relation with four parameters that provided excellent description of data patterns for a wide range of adsorption systems.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Weber Van Vliet isotherm analysis
#' @export
webervvanalysis <-function(Qe,Ce){
	x<- Qe
	y <- Ce

	fit81 <- lm(y~x)
	print(summary(fit81))
	rhs <- function(x,b0,b1,b2,b3){
	x^(b0*x^(b1) + b2)
	}
	fit82 <- nlsLM(y~ rhs(x,p2,p3,p4) , start=list( p2=1, p3=1, p4=1), trace = TRUE)
	print("Weber- Van Vliet Isotherm")
	print(summary(fit82))

	a <- (summary(fit82))
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
