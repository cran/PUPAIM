#' @title Elovich Kinetics
#' @description This equation assumes that the actual solid surfaces are energetically heterogenous and that neither desorption nor interactions betwen adsorbed species could subcutaneously affect the kinetics of adsorption at low surface coverage. (Mercado-Borayo, et. al, 2014)

#' @param t duration for the experiment
#' @param Qt the numerical value for the concentration at given time
#'
#' @return The regression analysis for the Elovich Kinetics
#'
#'@export
elovichkinetics <-  function(t,Qt){
  x <- log10(t)
  y <- Qt

  fit17 <- lm(y~x)
  print("Elovich Kinetics Analysis")
  print(summary(fit17))

  rhs <- function(x,b0,b1){
    1/(b0)*log10(b0*b1) + x
  }
   fit18 <- nlsLM(y~rhs(x,a,b), start=list(a= 0.0985913,b= 18.3412), trace = TRUE)
   print("Parameters")
   print(summary(fit18))
   a <- (summary(fit17))
   b <- a$coefficients[2]
   c <- a$coefficients[1]
   print(b)
   print(c)
   Qp <- function(qt){
     x <- Qt

     j<- (b*x) + c

     print(j)
   }
   b4 <- Qp(Qt)

   errors <- function(Qt,Qp){

     rmse<- rmse(Qt,b4)
     mae<- mae(Qt,b4)
     mse <- mse(Qt,b4)
     rae <-rae(Qt,b4)
     colnames(Qt) <- rownames(Qt) <-colnames(Qt)
     list("relative mean squared error" = rmse, "mean absolute error"=mae , "mean squared error"=mse ,  "relative absolute error"=rae)

   }

   a <- errors(Qt,Qp)
   print(a)


}

