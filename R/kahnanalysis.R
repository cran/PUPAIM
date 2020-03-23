#' @title Kahn Isotherm Analysis Non-Linear Form
#' @description Kahn isotherm model is a general model for adsorption of biadsorbate from pure dilute equations solutions
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value the absorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the linear regression and the parameters for the Kahn isotherm analysis
#' @examples kahnanalysis(c(0.5, 0.6, 0.7, 1.3, 3.9, 7.6, 16.5),c(4.5, 6.5, 9.3, 10.7, 11.6, 12.5, 13.7))
#' @author Jann Audrey B. Angeles
#' @author Raymond B. Diaz
#' @author C.C. Deocaris
#' @references A. R. Khan, R. Ataullah,A. Al-Haddad. (1997). Equilibrium Adsorption Studies of Some Aromatic Pollutants from Dilute Aqueous Solutions on Activated Carbonat Different Temperatures. JOURNAL OF COLLOID AND INTERFACE SCIENCE, 154-165.
#' @references Tosun, I. (2012). Ammonium Removal from Aqueous Solutions by Clinoptilolite: Determination of Isotherm and Thermodynamic Parameters and Comparison of Kinetics by the Double Exponential Model and Conventional Kinetic Models. International Journal of Environmental Research and Public Health , 970-984.
#' @export
kahnanalysis <- function(Ce,Qe) {
  x <- Ce
  y <- Qe
  data <- data.frame(Ce, Qe)
print("NLS2 Analysis for Kahn Isotherm")
mod240 <- (Qe ~ (Qmax*bk*Ce)/((1+bk*Ce)^ak))
start <- data.frame(Qmax = c(1, 10),
                    bk = c(1, 100),
                    ak = c(0, 1))
set.seed(511)
fit241 <- nls2(mod240, start = start,
           control = nls.control(maxiter =1000 , warnOnly = TRUE),
           algorithm = "plinear-random")
print(summary(fit241))
N <- nrow(na.omit(data))
error <- function(y){
  pv  <- (predict(fit241))
  rmse<- (rmse(y,predict(fit241)))
  mae <- (mae(y,predict(fit241)))
  mse <- (mse(y,predict(fit241)))
  rae <- (rae(y,predict(fit241)))
  PAIC <- AIC(fit241)
  PBIC <- BIC(fit241)
  SE <-(sqrt(sum(predict(fit241)-x)^2)/(N-2))
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
  plot(x,y,main = "Kahn Isotherm Plot", xlab = "Ce", ylab = "Qe")
  lines(smooth.spline(x, predict(fit241)), col = "black")
  rsqq <- lm(Qe~predict(fit241))
  print(summary(rsqq))
}
