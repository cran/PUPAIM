#' @title Kahn Isotherm
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Kahn isotherm analysis
#' @examples kahnanalysis(moringa$Ce,moringa$Qe)
#' @export
kahnanalysis <-function(Ce,Qe){
	x <- Ce
	y <- Ce/Qe
	fit53<- lm(y~x)
	print("Kahn Analysis")
	print(summary(fit53))

	rhs <- function(x, b0,b1){
	b0 + b1*b0*x
	}
	fit54 <- nlsLM(y~rhs(x,ak,bk), start=list(ak=1, bk=1), trace=TRUE)
	print("Parameters")
	print(summary(fit54))

	a <- (summary(fit53))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit54))
	e <- d$coefficients[1]
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
	deltaG(25,e)
	}
