#' @title Hill Deboer Isotherm Analysis
#' @description Describes an incident where there is mobile adsorption as well as lateral interaction among adsorbed molecules.
#' @param theta the numerical value representing the fractional coverage
#' @param Ce the numerical value for the equilibrium capacity
#'
#' @return the linear regression and the parameters for the Hill Deboer isotherm analysis
#' @examples hilldeboeranalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
hilldeboeranalysis <-function(theta, Ce){
	x <- theta/24.46519
  y <- -log10(Ce)
 fit43 <- lm(y~x)
print("Hill-Deboer Analysis")
print(summary(fit43))
 rhs <- function(x, b0, b1){
   -log10(b0) - b1*x
 }
  fit44 <- nlsLM(y ~ rhs(x, k1,k2), start=list(k1=1, k2=1), trace=TRUE)
  print(summary(fit44))
  a <- (summary(fit43))
  b <- a$coefficients[2]
  c <- a$coefficients[1]
  d <-(summary(fit44))
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
