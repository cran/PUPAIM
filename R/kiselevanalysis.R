#' @title Kiselev Isotherm Analysis
#'
#' @param Kh the numerical value representing Kiselev analysis
#' @param Ce the numerical value for the equilibrium capacity
#'
#' @return the linear regression and the parameters for the Kiselev isotherm analysis
#' @examples kiselevanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
kiselevanalysis <-function(Kh,Ce){

	x<- Kh

	y<-1/Ce

	fit56 <- lm(y~x)
	print(summary(fit56))

  rhs<-function(x, b1, b2){

	1/b1 + (b2*x)

}

	fit57<-nlsLM(y~rhs(x, K1, Ki), start=list (K1= 244.983, Ki=3.7047e+11), trace=TRUE)

	print(summary(fit57))
	a <- (summary(fit56))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit57))
	e <- d$coefficients[1]
	print(b)
	print(c)
	Qp <- function(Kh){
	  x <- Kh

	  j<- (b*x) + c

	  print(j)
	}
	b4 <- Qp(Kh)

	errors <- function(Kh,Qp){

	  rmse<- rmse(Kh,b4)
	  mae<- mae(Kh,b4)
	  mse <- mse(Kh,b4)
	  rae <-rae(Kh,b4)
	  colnames(Kh) <- rownames(Kh) <-colnames(Kh)
	  list("Predicted values", "relative mean squared error" = rmse, "mean absolute error"=mae , "mean squared error"=mse ,  "relative absolute error"=rae)}
	errors(Kh,Qp)

	deltaG(25,e)
}

