#' @title Langmuir Isotherm Analysis
#' @description It was primarily designed to describe gas-solid phase adsorption is also used to quantify and contrast the adsorptive capacity of various adsorbents [12]. Langmuir isotherm accounts for the surface coverage by balancing the relative rates of adsorption and desorption (dynamic equilibrium).
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom Metrics "rmse"
#' @importFrom Metrics "mae"
#' @importFrom Metrics "rae"
#' @importFrom Metrics "mse"
#' @return the regression analysis for the langmuir isotherm analysis
#' @examples langmuiranalysis(moringa$Ce,moringa$Qe)
#' @export
langmuiranalysis <- function (Ce, Qe) {

    x1 <- 1/Ce
    y1 <- 1/Qe


	rhs <- function(x,b0,b1){
	  (1/b0*b1)*x1 + (1/b0)

	}
	fit62<- nlsLM(y1~rhs(x1,Qm,K),start=list(Qm= 1,K=1))
	print(summary(fit62))

	fit63 <- lm(y1 ~ x1)
	print("Langmuir Isotherm Analysis")
	print(summary(fit63))

	a <- (summary(fit63))
b <- a$coefficients[2]
c <- a$coefficients[1]
d <- (summary(fit62))
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
  list("relative mean squared error" = rmse, "mean absolute error"=mae , "mean squared error"=mse ,  "relative absolute error"=rae)

}

 errors(Ce,Qp)
deltaG(25,e)





}



