#' @title Fritz Schlunder Isotherm Analysis
#' @description An empirical equation which can fit a wide range of experimental results because of the large number of coefficients in the isotherm.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the linear regression and the parameters for the Fritz Schlunder isotherm analysis
#' @examples fritzanalysis(moringa$Ce,moringa$Qe)
#' @export
#'
fritzanalysis <-function(Ce,Qe){

	x <- Ce

	y <- Qe

	fit28 <- lm(y~x)

	print("Fritz-Schlunder Analysis")

	print(summary(fit28))

	rhs <- function(x, b0, b1, b2){

    	(x*b1)/(1+b2*x^b0)

}

	fit29 <- nlsLM(y~ rhs(x, Qm, Kfs, Mfs), start=list (Qm=1, Kfs=1, Mfs=1), trace=TRUE)

	print("Fritz-Schlunder Parameters")

	print(summary(fit29))

	a <- (summary(fit29))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit29))
	e <- d$coefficients[2]
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
