#' @title Fowler Guggenheim Isotherm
#' @description This isotherm equation which takes into consideration the lateral interaction of the adsorbed molecules.
#' @param theta the numerical value reperesenting the fractional coverage
#' @param Ce the numerical value fot the equilibrium capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the linear regression and the parameters for the Fowler Guggenheim isotherm analysis
#' @examples fowlerganalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
fowlerganalysis <-function(theta, Ce){
  x<- theta/24.46529
  y <- -log10(Ce)
  fit23 <- lm(y~x)
  print("Fowler Guggenheim Analysis")
  print(summary(fit23))


  rhs<- function(x,b0,b1){
    -log10(b0)+(2*b1*x)
  }
  fit24 <- nlsLM(y~ rhs(x,K,w), start=list(K=1, w=1), trace=TRUE)
  print("Fowler Guggenheim Parameters")
  print(summary(fit24))

  a <- (summary(fit23))
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
