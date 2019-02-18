#' @title Marckzewski Jaroniec Isotherm Analysis
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Marckzewski Jaroniec isotherm analysis
#'
#' @export
mjanalysis <-function(Ce,Qe){

	x<- 1/Ce^2

	y<-Qe

	fit67 <- lm(y~x)
	print("MJ Analysis")
	print(summary(fit67))

	rhs<-function(x, b1, b2,b3){
	  ((b1)/(1+(b1*x)))^(b2/b3)

}

	fit68<-nlsLM(y~rhs(x, Kmj, nmj,Mmj),start=list(Kmj=2.31123e+06, nmj=6325.36,Mmj=3319.24), trace=TRUE)
  print("Parameters")
	print(summary(fit68))
	a <- (summary(fit67))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit68))
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
