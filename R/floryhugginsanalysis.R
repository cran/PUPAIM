#' @title Flory Huggins Isotherm Analysis
#'
#' @param theta the numeric value for the fractional coverage
#' @param Co the numeric value for the initial concentration
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "nls"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the regression analysis for the Flory Huggins Isotherm
#' @export
fhanalysis <- function(theta,Co){

  	x2 <- log10(1-theta)

	y2 <- log10(theta/Co)

	fit20 <- lm(y2~x2)

	print("Flory-Huggins Analysis")

	print(summary(fit20))

	rhs <- function(x2, b0, b1){

(-log10(b0)+b1*x2)

	}

	fit21 <- nls(y2~ rhs(x2, Kfh, n), start=list (Kfh=0.5, n=1))

	print("Flory-Huggins Parameters")

	print(summary(fit21))
	a <- (summary(fit20))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit21))
	e <- d$coefficients[1]
	print(b)
	print(c)
	Qp <- function(theta){
	  x <- theta

	  j<- (b*x) + c

	  print(j)
	}
	b4 <- Qp(theta)

	errors <- function(theta,Qp){

	  rmse<- rmse(theta,b4)
	  mae<- mae(theta,b4)
	  mse <- mse(theta,b4)
	  rae <-rae(theta,b4)
	 colnames(theta) <- rownames(theta) <-colnames(theta)
	  list("Predicted values", "relative mean squared error" = rmse, "mean absolute error"=mae , "mean squared error"=mse ,  "relative absolute error"=rae)}
	errors(theta,Qp)
	deltaG(25,e)
}
