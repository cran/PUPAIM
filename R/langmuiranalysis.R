#' @title Langmuir Isotherm Analysis Non-Linear Form
#' @description The Langmuir adsorption isotherm is used to describe the equilibrium between adsorbate and adsorbent system, where the adsorbate adsorption is limited to one molecular layer at or before a relative pressure of unity is reached.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the linear regression and the parameters for the Langmuir isotherm
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples langmuiranalysis (Ce, Qe)
#' @author Mark Lester Galicia
#' @author C.C. Deocaris
#' @references Langmuir, I. (1997). The constitution and fundamental properties of solids and liquids II. Liquids. Journal of American Chemistry Society, 1848-1906.
#' @references Liu, L., Luo, X.-B., LinDing, & Sheng-LianLuo. (2019). Application of Nanotechnology in the Removal of Heavy Metal From Water. Nanomaterials for the Removal of Pollutants and Resource Reutilization, 83-147.
#' @export
langmuiranalysis <- function(Ce, Qe){
  x<-Ce
  y<-Qe
  data <- data.frame(x, y)
  mod247<-(Qe~(Qmax*b*Ce)/(1+(b*Ce)))
  set.seed(511)
  start<-data.frame(Qmax = c(0,1000), b = c(0,1000))
  suppressWarnings(fit248<-nls2(mod247, start = start,
                               control = nls.control(maxiter=50, warnOnly = TRUE),
                               algorithm="port"))
  print("Langmuir Isotherm Analysis")
  print(summary(fit248))
  N<- nrow(na.omit(data))
  error <- function(y){
    pv  <- (predict(fit248))
    rmse<- (rmse(y,predict(fit248)))
    mae <- (mae(y,predict(fit248)))
    mse <- (mse(y,predict(fit248)))
    rae <- (rae(y,predict(fit248)))
    PAIC <- AIC(fit248)
    PBIC <- BIC(fit248)
    SE <-(sqrt(sum(predict(fit248)-x)^2)/(N-2))
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
  plot(x,y, main = "Langmuir Non-Linear Isotherm Model", xlab = "Ce",
       ylab = "Qe")
  lines(data$x, predict(fit248), col="black", lwd =1)
  rsqq <- lm(Qe~predict(fit248))
  print(summary(rsqq))
}
