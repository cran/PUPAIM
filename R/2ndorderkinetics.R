#' @title Pseudo-2nd Order Kinetics
#' @description The pseudo-second order model describes the adsorption reaction rate with depndent energetically hetrogenousites on the adsorbent. (Mercado-Borayo, et. al., 2014)
#' @param t duration of the experiment
#' @param Ce the numerical value for the equilibrium capacity
#' @return the regression analysis for the second order kinetics
#' @importFrom graphics "plot"
#' @importFrom graphics "abline"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @importFrom minpack.lm "nlsLM"
#' @examples secondorder(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
secondorder <- function(t,Ce){
  x <- t
  y <- 1/Ce
  fit3 <- lm(y~x)
  print("Second order Kinetics")
  print(summary(fit3))
  rhs <- function(x,b0,b1){
    1/(b1) +1/(b0*x^2)
  }
    fit4 <- nlsLM(y~rhs(x,K,Qe), start= list(K=1, Qe=1))
print("Parameters")
        print(summary(fit4))
        plot(y~x)
abline(fit3)
a <- (summary(fit3))
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
        }

