#' @title Henry Isotherm
#' @description It describes an appropriate fit to the adsorption of adsorbate at relatively low concentrations such that all adsorbate molecules are secluded from their nearest neighbours.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Henry isotherm analysis
#' @examples henryanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
henryanalysis <-function(Ce, Qe){
	x <- Ce
	y <- Qe
	fit38 <- lm(y ~ x)
	print("Henry Isotherm Analysis")
	print(summary(fit38))

	rhs <- function(x,b0){
	  b0*x
	}
	 fit39 <- nlsLM(y~ rhs(x, K), start=list(K=1), trace=TRUE)
	 print("Parameters")
	 print(summary(fit39))


	 a <- (summary(fit38))
	 b <- a$coefficients[2]
	 c <- a$coefficients[1]
	 d <- (summary(fit39))
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
