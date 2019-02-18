#' @title Sips Isotherm Analysis
#'
#' @param Ce the numerical value for the equilbrium concentration
#' @param Qe the numerical value for the adsorbed concentration
#'
#' @return the Linear model for the Sips isotherm analysis
#' @examples sipsanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#'
#' @export
sipsanalysis <-function(Qe,Ce){
	x <- -log(1/Qe)
	y <- log(Ce)
	fit76 <- lm(y~x)
	print("Sips Isotherm Analysis")
	print(summary(fit76))





	a <- (summary(fit76))
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
