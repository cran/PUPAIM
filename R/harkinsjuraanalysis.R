#' @title HarkinsJura Isotherm
#' @description It assumes the possibility of multilayer adsorption on the surface of absorbents having heterogeneous pore distribution.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#' @return the linear regression and the parameters for the Harkins Jura isotherm analysis
#' @examples harkinsjuraanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
harkinsjuraanalysis <-function(Ce,Qe){
	x<-log(Ce)
	y<- (1/(Qe)^2)
	fit35<-lm(y~x)
	print("Harkins-Jura Isotherm Analysis")
	print(summary(fit35))

	rhs <- function(x,b1){
	  (1/b1)*x
	}
	fit36 <- nlsLM(y~rhs(x,A), start=list(A= 1), trace=TRUE)
	print("Parameters")
	print(summary(fit36))
	a <- (summary(fit35))
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
