#' @title Redlich Peterson Isotherm Analysis
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Redlich Peterson isotherm analysis
#' @examples redlichpanalysis(moringa$Ce,moringa$Qe)
#' @export
redlichpanalysis <-function(Ce, Qe) {

	x <- log10(Ce)
	y <- log10 (Ce/Qe)

	fit73 <- lm(y ~ x)
	print ("Redlich Peterson Analysis")
	print (summary(fit73))

	rhs <- function(x, b0){
	  b0*x
		}
			fit74 <- nls(y~ rhs(x, g), start= list(g=1.0))
	print("Redlich Peterson Parameters")
	print(summary(fit74))


	a <- (summary(fit73))
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
