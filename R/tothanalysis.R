#' @title Toth Isotherm
#' @description Another empirical modification of the Langmuir equation with the aim of reducing the error between experimental data and predicted value of equilibrium data.
#' @param Ce the numerical value for the equilibrium capacity
#' @param theta the numerical value representing the fractional coverage
#'
#' @return the linear regression and the parameters for the Toth isotherm analysis
#' @export
tothanalysis <-function(Ce,theta){

	x <- Ce

	y <- theta

	fit78 <- lm(y~x)
	print(summary(fit78))
	rhs<-function(x, b0, b1,b2){

	(b0*x)/(1+(b1*x))^b2

}

	fit79 <- nlsLM(y~rhs(x, Ke, Kl,n), start=list(Ke=1 , Kl=1,n=2), trace=TRUE)
	print("Toth Isotherm")
	print(summary(fit79))


	a <- (summary(fit78))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit79))
	e <- d$coefficient[2]
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
