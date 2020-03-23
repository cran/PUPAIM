#' @title Langmuir Isotherm Analysis Linear Form
#' @description The Langmuir Linear Equation I adsorption isotherm is used to describe the equilibrium between adsorbate and adsorbent system, where the adsorbate adsorption is limited to one molecular layer at or before a relative pressure of unity is reached.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "plot" "abline"

#' @importFrom minpack.lm "nlsLM"
#' @return the linear regression and the parameters for the Langmuir isotherm
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples langmuir.LM1(Ce, Qe)
#' @author Mark Lester Galicia
#' @author C.C. Deocaris
#' @references Langmuir, I. (1917). The constitution and fundamental properties of solids and liquids II. Liquids. Journal of American Chemistry Society, 1848-1906.
#' @references Liu, L., Luo, X.-B., LinDing, & Sheng-LianLuo. (2019). Application of Nanotechnology in the Removal of Heavy Metal From Water. Nanomaterials for the Removal of Pollutants and Resource Reutilization, 83-147.
#' @export
langmuir.LM1 <- function (Ce, Qe)
{
  x <- 1/Ce
  y <- 1/Qe
  dat <- data.frame(x,y)
  N<- nrow(na.omit(dat))
  rhs <- function(x, b, Qmax) {
    (1/(b * Qmax)) * x + (1/Qmax)
  }
  fit251 <- lm(y ~ x)
  print("Langmuir Isotherm Analysis")
  print(summary(fit251))
  a <- (summary(fit251))
  b <- a$coefficients[2]
  c <- a$coefficients[1]
  param.ads <- function(y){
    qmax <- 1/a$coefficients[1]
    b <- (a$coefficients[2]/qmax)
    colnames(y) <- rownames(y) <- colnames(y)
    print("Langmuir Parameters")
    list("Qmax" = qmax,
         "b" = b)
  }
  w <- param.ads(y)
  print(w)
  error <- function(y){
    pv  <- (predict(fit251))
    rmse<- (rmse(y,predict(fit251)))
    mae <- (mae(y,predict(fit251)))
    mse <- (mse(y,predict(fit251)))
    rae <- (rae(y,predict(fit251)))
    PAIC <- AIC(fit251)
    PBIC <- BIC(fit251)
    SE <-(sqrt(sum(predict(fit251)-x)^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Predicted Values"           = pv,
         "Relative Mean Square Error" = rmse,
         "Mean Absolute Error"        = mae,
         "Mean Squared Error"         = mse,
         "Relative Absolute Error"    = rae,
         "AIC"                        = PAIC,
         "BIC"                        = PBIC,
         "Standard Error Estimate"    = SE)
  }
  e <- error(y)
  print(e)
  plot(x, y, main = "Langmuir Isotherm Linear Model 1", xlab = "1 / Ce",
       ylab = "1 / Qe")
  abline(fit251, col="black")
}
