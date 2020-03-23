#' @title Temkin Isotherm Analysis Linear Form
#' @description takes into account the effects of indirect adsorbate/ adsorbate interaction on the adsorption process
#' @param Ce wherein theta is  is the numerical value for the equilibrium capacity
#' @param Qe is Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Flory-Huggins isotherm
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples temkinanalysis(Ce, Qe)
#' @author Lance E. Abadier
#' @author C.C. Deocaris
#' @references Newton, A., & Donbebe. (2017, September 5). Modelling and Interpretation of Adsorption Isotherms. retrieved from: https://www.hindawi.com/journals/jchem/2017/3039817
#' @export
temkinanalysis <- function(Ce, Qe){
  x <- Ce
  y <- Qe
  data<- data.frame(Ce, Qe)
  n<- nrow(na.omit(data))
  fit262 <- Qe ~ (B*log10(A*Ce))
  print("Temkin Analysis")
  start <- data.frame(A= c(0, 1000), B= c(0, 100))
  set.seed(511)
  fit263 <- nls2(fit262, start = start, control = nls.control(maxiter = 100, warnOnly = FALSE), algorithm = "port")
  print(summary(fit263))
  error <- function(y){
    pv  <- (predict(fit263))
    rmse<- (rmse(y,predict(fit263)))
    mae <- (mae(y,predict(fit263)))
    mse <- (mse(y,predict(fit263)))
    rae <- (rae(y,predict(fit263)))
    PAIC <- AIC(fit263)
    PBIC <- BIC(fit263)
    SE <-(sqrt(sum(predict(fit263)-x)^2)/(n-2))
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
  plot(x, y, main = "Temkin Isotherm Non-Linear Plot", xlab = "Ce", ylab = "Qe")
  lines(x, predict(fit263), col = "black")
  rsqq <- lm(y~predict(fit263))
  print(summary(rsqq))
}
