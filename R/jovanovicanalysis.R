#' @title Jovanovic Isotherm
#' @description It is predicated on the assumptions contained in the Langmuir model, but in addition the possibility of some mechanical contacts between the adsorbate and adsorbent.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Jovanovic isotherm analysis
#' @examples jovanovicanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
jovanovicanalysis <-function(Ce,Qe){
	x5<-Ce
	y5<-log10(Qe)
	fit50 <-lm(y5~x5)
	print("Jovanovic Isotherm Analysis")
	print(summary(fit50))

	 rhs<- function(x5, b0,b1){
	-b0*x5
}
		fit51<- nls(y5 ~ rhs(x5, Kj), start=list(Kj=1))
	print("Parameters")
		print(summary(fit51))
	a <- (summary(fit50))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit51))
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
