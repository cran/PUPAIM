#' @title Pseudo-1st Order Kinetics
#' @description A first-order rate equation which is believed to be the earliest model that was presented by Lagergen (1898) to describe the kinetic process of liqud-solid phase adsorption of oxalic acid and malonic acid onto charcoal pertains to the adsorption rate based on the adsorption capacity.
#' @param t duration of the experiment
#' @param Ce the numerical value for the equilibrium capacity
#'
#' @return the regression analysis for the first order kinetics
#' @examples firstorder(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
firstorder <- function(t,Ce){
  x1 <- t
  y1 <- log10(Ce)
  fit <- lm(y1~x1)
  print("First order Kinetics")
  print(summary(fit))

  rhs <- function(x,b0,b1){
    -(b0*x) + log10(b1)

  }

   fit1 <- nls(y1~rhs(x1, K,Qe), start= list(K=1, Qe=1))
   print("Parameters")
   print(summary(fit1))
   plot(y1~x1)
   abline(fit)


   a <- (summary(fit))
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
