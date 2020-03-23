#' @title Jossens Isotherm Analysis Non-Linear Form
#' @description The Jossens isotherm model predicts a simple equation based on the energy distribution of adsorbate-adsorbent interactions at adsorption sites. This model assumes that the adsorbent has heterogeneous surface with respect to the interactions it has with the adsorbate.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Jossens isotherm
#' @examples jossensanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @author Stephanie Mae L. Manuel
#' @author Marie Aileen M. Allauigan
#' @author C.C. Deocaris
#' @references Ayawei, N., Ebelegi, A., & Wankasi, D. (2017, September 5). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 6. doi:10.1155/2017/3039817
#' @export
jossensanalysis <- function(Qe, Ce){
x <- Qe
y <- Ce
data <- data.frame(x, y)
print("NLS2 Analysis for Jossens Isotherm")
mod236 <- (Ce ~ (Qe/b0)*(exp(b1*(Qe^b2))))
start <- data.frame(b0 = c(0, 500), b1 = c(0.0001, 100), b2 = c(0, 1))
set.seed(511)
fit237 <- nls2(mod236, start = start, control = nls.control(maxiter =100 , warnOnly = TRUE), algorithm = c("plinear-random"))
print(summary(fit237))
N <- nrow(na.omit(data))
error <- function(y){
  pv  <- (predict(fit237))
  rmse<- (rmse(y,predict(fit237)))
  mae <- (mae(y,predict(fit237)))
  mse <- (mse(y,predict(fit237)))
  rae <- (rae(y,predict(fit237)))
  PAIC <- AIC(fit237)
  PBIC <- BIC(fit237)
  SE <-(sqrt(sum(predict(fit237)-x)^2)/(N-2))
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
  plot(x, y, main ="Jossens Isotherm Plot", xlab= "Qe", ylab= "Ce")
  lines(x, predict(fit237), col = "black")
  rsqq <- lm(y~predict(fit237))
  print(summary(rsqq))
}
