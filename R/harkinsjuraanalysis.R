#' @title HarkinsJura Isotherm Analysis Non-Linear Form
#' @description It assumes the possibility of multilayer adsorption on the surface of absorbents having heterogenous pore distribution (Ayawei, et al., 2017),(Gupta, et al., 2012)
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the HarkinsJura isotherm
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples harkinsjuraanalysis(Ce, Qe)
#' @author Raymond James L. Adame
#' @author C.C. Deocaris
#' @references Ayawei, N., Ebelegi, A. N., & Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 1-11. doi: 10.1155/2017/3039817
#' @references Gupta, V. K., Pathania, D., Agarwal, S., & Sharma, S. (2012). Removal of Cr(VI) onto Ficus carica biosorbent from water. Environmental Science and Pollution Research, 20(4), 2632-2644. doi:10.1007/s11356-012-1176-6
#' @export
harkinsjuraanalysis <- function(Ce, Qe){
x <- Ce
y <- Qe
data <- data.frame(x, y)
print("Harkins-Jura Isotherm Analysis")
set.seed(511)
mod227 <- (y ~ (A/(b-log(x)))^1/2)
start <- data.frame(A = c(0, 100), b = c(0, 100))
suppressWarnings(fit228 <- nls2(mod227, start = start, control = nls.control(maxiter =50 , warnOnly = TRUE), algorithm = "port"))
fit228
print(summary(fit228))
N <- nrow(na.omit(data))
error <- function(y){
  pv  <- (predict(fit228))
  rmse<- (rmse(y,predict(fit228)))
  mae <- (mae(y,predict(fit228)))
  mse <- (mse(y,predict(fit228)))
  rae <- (rae(y,predict(fit228)))
  PAIC <- AIC(fit228)
  PBIC <- BIC(fit228)
  SE <-(sqrt(sum(predict(fit228)-x)^2)/(N-2))
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
plot(x, y, main = "Harkins-Jura Isotherm", xlab = "Ce", ylab = "Qe")
lines (x, predict(fit228), col = "black")
rsqq <- lm(Qe~predict(fit228))
print(summary(rsqq))
}
