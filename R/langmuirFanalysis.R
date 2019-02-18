#' @title Langmuir Freundlich Isotherm Analysis
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the analysis for the Langmuir-Freundlich Isotherm
#' @examples langmuirFanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
langmuirFanalysis <- function(Ce, Qe){
	x<- 1/(Ce)^2
	y<- Qe

	fit64  <-lm(y~x)

	print(summary(fit64))
	rhs <- function(x,b0,b1,b2){
		b2*((1+b0)^b1*x)/(1+(1+b0*x)^b1)
	}
    fit65<-nlsLM(y~rhs(x,f,g,h), start=list(f=1,g=1,h=1), trace=TRUE)

	  print(summary(fit65))
	  a <- (summary(fit64))
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
