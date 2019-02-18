#' @title Raudke Prausniiz Isotherm Analysis
#' @description  It has several important properties which makes it more preferred in most adsorption systems at low adsorbate concentration.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Raudke Prausniiz isotherm analysis
#' @examples raudkepanalysis(moringa$Ce,moringa$Qe)
#' @export
raudkepanalysis <-function(Ce, Qe) {

	x1 <- Ce^2
	y1 <- Ce/Qe

	fit70<- lm(y1~x1)
	print ("Raudke-Prausniiz Analysis")
	print(summary(fit70))



	rhs<- function(x, b0, b1){
	  (1+(b0*x)^(b1))
	}
	fit71 <- nlsLM(y1~ rhs(x1, Krp, MRP), start= list(Krp=0.85 , MRP=1), trace=TRUE)
	print("Raudke-Prausniiz Parameters")
	print(summary(fit71))
	a <- (summary(fit70))
	b <- a$coefficients[2]
	c <- a$coefficients[1]
	d <- (summary(fit71))
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
