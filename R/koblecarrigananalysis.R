#' @title Koble Carrigan Isotherm Analysis
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the linear regression and the parameters for the Koble Carrigan isotherm analysis
#' @examples koblecarrigananalysis(moringa$Ce,moringa$Qe)
#' @export
koblecarrigananalysis <-function(Ce, Qe){
 x <- 1/Ce
 y <- 1/Qe

 fit59 <- lm(y ~x)
 print("Koble Carrigan Isotherm")
 print(summary(fit59))

 rhs <- function( x, b0,b1){
   1/b0*x + (b1/b0)
 }
 fit60 <- nlsLM(y~rhs(x, Ak, Bk) , start=list(Ak=1, Bk=1))
 print("Parameters")
 print(summary(fit60))
 a <- (summary(fit59))
 b <- a$coefficients[2]
 c <- a$coefficients[1]
 d <-(summary(fit60))
 e<- d$coefficients[1]
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
   list("predicted values", "relative mean squared error" = rmse, "mean absolute error"=mae , "mean squared error"=mse ,  "relative absolute error"=rae)

 }

 a <- errors(Ce,Qp)
 print(a)
deltaG(25,e)



 }
