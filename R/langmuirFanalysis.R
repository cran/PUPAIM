#' @title Langmuir-Freundlich Isotherm Analysis Non-Linear Form
#' @description Langmuir-Freundlich Isotherm Analysis describes the distribution of adsorption energy onto heterogeneous surfaces of the adsorbent. At low concentrations of adsorbate, this model becomes Freundlich isotherm, and then at high concentration of adsorbate, it becomes the Langmuir Isotherm. The parameters of this concentration can be obtained using the non-linear regression.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the linear regression and the parameters for the Langmuir-Freundlich isotherm
#' @examples Ce <- c(0.06649, 0.21948, 0.38188, 0.56311, 0.77729, 0.98794, 1.25390, 1.72698)
#' @examples Qe <- c(0.05192, 0.07174, 0.08680, 0.09433, 0.08839, 0.08363, 0.09711, 0.10741)
#' @examples langmuirFanalysis(Ce, Qe)
#' @author Cabugnason, Jay Anne M.
#' @author Pogado, Precious Grace
#' @author C.C. Deocaris
#' @references Jeppu G.P and Clement T.P. (2012, March 15) A modified Langmuir-Freundlich isotherm model for simulating pH-dependent adsorption effects. Retrieved from: www.researchgate.net/publication/221762917_A_modified_Langmuir-Freundlich_isotherm_model_for_simulating_pH-dependent_adsorption_effects
#' @references Umpleby R.J., et. al (2001, August 23) Characterization of Molecularly Imprinted Polymers with the Langmuir-Freundlich Isotherm. Retrieved fromhttps://pubs.acs.org/doi/pdf/10.1021/ac0105686
#' @export
langmuirFanalysis <- function(Ce, Qe) {
  x <- Ce
  y <- Qe
  data <- data.frame(x, y)
  print("Langmuir-Freundlich Isotherm Analysis")
  mod249 <- (y ~(Qmax*((klf*x)^n))/(1+((klf*x)^n)))
  start <- data.frame(Qmax = c(0, 1000),
                      klf = c(0, 1000),
                      n=c(0,1))
  set.seed(511)
  fit250 <- nls2(mod249,
                start = start,
                control = nls.control(maxiter = 100 , warnOnly = FALSE),
                algorithm = "plinear-random",trace=FALSE)
 print(summary(fit250))
 N <- nrow(na.omit(data))
 error <- function(y){
   pv  <- (predict(fit250))
   rmse<- (rmse(y,predict(fit250)))
   mae <- (mae(y,predict(fit250)))
   mse <- (mse(y,predict(fit250)))
   rae <- (rae(y,predict(fit250)))
   PAIC <- AIC(fit250)
   PBIC <- BIC(fit250)
   SE <-(sqrt(sum(predict(fit250)-x)^2)/(N-2))
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
plot(x, y, main = "Langmuir-Freundlich Isotherm Plot", xlab = "Ce", ylab = "Qe")
lines(x, predict(fit250), col = "black")
rsqq <- lm(Qe~predict(fit250))
print(summary(rsqq))
}
