#' @title Baudu Isotherm Analysis
#' @description Reduced form of Langmuir Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'@return the linear regression and the parameters for the Baudu isotherm
#'@examples bauduanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#'
#' @export
#'
bauduanalysis <- function(Ce, Qe) {

	x1 <- 1/Ce
	y1 <- Qe

	fit8 <- lm(y1 ~ x1)
	print ("Baudu Analysis")
	print (summary(fit8))

	rhs <- function(x1, b0, b1){
  (b0*b1*x1)/(1+(b0*x1))
	  }
			fit9 <- nlsLM(y1~ rhs(x1, qm, B), start= list( qm=28.7461, B=39.9559), trace=TRUE)

	print("Baudu Parameters")
	print(summary(fit9))
	a <- (summary(fit8))
	b <- a$coefficients[2]
	c <- a$coefficients[1]

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
	v <- errors(Ce,Qp)
print(v)
}

