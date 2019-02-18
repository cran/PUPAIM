#' @title Freundlich Isotherm Analysis
#' @description Gives an expression which defines the surface heterogeneity and the exponential distribution of active sites and their energies
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the regression analysis for the freundlich isotherm analysis
#' @examples freundlichanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
#'
#'
freundlichanalysis <- function (Ce, Qe)
{
    x2 <- log10(Ce)
    y2 <- log10(Qe)
rhs <- function(x,b0,b1){
log10(b0)+(1/b1)*x2
}
 fit26 <- nlsLM(y2~rhs(x2,Kf,n),start=list(Kf=1,n=1))
	print(summary(fit26))
	fit27 <- lm(y2 ~ x2)

	print("Freundlich Isotherm Analysis")
	print(summary(fit27))

	a <- (summary(fit27))
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
	  list("relative mean squared error" = rmse, "mean absolute error"=mae , "mean squared error"=mse ,  "relative absolute error"=rae)

	}

	a <- errors(Ce,Qp)
	print(a)
deltaG(25,b)

}
