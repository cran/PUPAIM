#' @title Five Parameter Isotherm Analysis
#' @description A five-parameter empirical model that is capable of simulating the model variations more precisely for application over a wide range of equilibrium data.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#'
#' @return the linear regression and the parameters for the five parameter isotherm analysis
#' @examples fiveparamanalysis(moringa$Ce,moringa$Qe)
#'@export
fiveparamanalysis <- function(Ce,Qe){
  x <- Ce^2
  y <- Qe

  fit5 <-lm(y~x)
  print("Five parameter analysis")
  print(summary(fit5))


  rhs <-function(x,b0,b1,b2,b3){
    (b0*x^b3)/(1+b1*x^b2)
  }

  fit6 <-nlsLM(y~rhs(x,k1,alpha,k2,beta),start = list(k1= 1,alpha=1,k2=1,beta=1),trace = TRUE)
  print(summary(fit6))
  a <- (summary(fit5))
  b <- a$coefficients[2]
  c <- a$coefficients[1]
  d <- (summary(fit6))
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
