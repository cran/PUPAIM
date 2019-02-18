#' @title Jossens Isotherm
#' @description It predicts a simple equation based on the energy distribution of adsorbate-adsorbent interactions at adsorption sites.This model assumes that the adsorbent has heterogeneous surface with respect to the interactions it has with the adsorbate.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Jossens isotherm analysis
#' @examples jossensanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
jossensanalysis <-function(Ce,Qe){
x<- Qe
y <-log10(Ce/Qe)

fit47 <- lm(y~x)

print("Jossen's Isotherm")
print(summary(fit47))


rhs <- function(x, b0,b1){
b0*x^b1
}
fit48 <- nlsLM(y~rhs(x, f, p), start=list(f=1, p=1))
print("Jossen's Parameters")
print(summary(fit48))
a <- (summary(fit47))
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
