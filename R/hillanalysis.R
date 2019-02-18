#' @title Hill Isotherm

#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Hill isotherm analysis
#' @examples hillanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
hillanalysis <-function(Ce, Qe){
	x <- log(Ce)
	y <- log(Qe/.5-Qe)
	fit41 <- lm(y~x)
	print("Hill Isotherm")
	print(summary(fit41))

	rhs<- function(x,b0,b1){
	x*b0 - log(b1)
		}

	fit42<- nlsLM(y~ rhs(x,nh,kD), start=list(nh=1, kD=1))
	print("Hill Parameters")
	print(summary(fit42))


	a <- (summary(fit41))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit42))
	e <- d$coefficients[2]
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
	  list("Predicted values", "relative mean squared error" = rmse, "mean absolute error"=mae , "mean squared error"=mse ,  "relative absolute error"=rae)
	  }
	v <- errors(Ce,Qp)
	print(v)
	print(deltaG(25,e))
	}
