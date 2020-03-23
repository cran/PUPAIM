#' @title Dubinin-Radushkevich Isotherm Analysis Non-Linear Form
#' @description It is an empirical generally applied to express adsorption mechanism with Gaussian energy distribution onto heterogeneous surfaces (Tsamo, et,al, 2019). In the non-linear form, qe =Xm*exp(-K*Ce^2). The said model has often fitted high solute activities and the intermediate range of concentrations as well. One unique feature of the said model is that it is temperature-dependent (Dada, et, al, 2012).
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "plot"
#' @importFrom graphics "abline"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the linear regression, errors, and the parameter values (Xm, and K) for the Dubinin-Radushkevich Isotherm Model
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples dubininradanalysis(Ce, Qe)
#' @author Morales, Dale Jonathan, M.
#' @author C.C. Deocaris
#' @references One-, Two-, and Three-Parameter Isotherms, Kinetics, and Thermodynamic Evaluation of Co(II) Removal from Aqueous Solution Using Dead Neem Leaves, Tsamo, et,al. (2019), Volume 2019, https://doi.org/10.1155/2019/6452672
#' @references Langmuir, Freundlich, Temkin, and Dubinin-Radushkevich Isotherms Studies of Equilibrium Sorption of Zn^+2 Unto Phosphoric Acid Modified Rise Husk, Dada, et, al, (2012), Journal of Applied Chemistry, Volume 3, Issue 1, pp 38-45, www.iosrjournals.org
#' @export
dubininradanalysis <- function(Ce, Qe){
  x <- Ce
  y <- Qe
  data <- data.frame(x, y)
  mod205 <- Qe ~ Xm*exp(-K*Ce^2)
  start <- data.frame(Xm = c(0, 10), K = c(0, 10))
  set.seed(511)
  suppressWarnings(fit206<- nls2(mod205, start = start, control = nls.control(maxiter = 100, warnOnly = TRUE), algorithm = "plinear-random"))
  print("Dubinin-Radushkevich Isotherm")
  print(summary(fit206))
  N <- nrow(na.omit(data))
  error <- function(y){
    pv  <- (predict(fit206))
    rmse<- (rmse(y,predict(fit206)))
    mae <- (mae(y,predict(fit206)))
    mse <- (mse(y,predict(fit206)))
    rae <- (rae(y,predict(fit206)))
    PAIC <- AIC(fit206)
    PBIC <- BIC(fit206)
    SE <-(sqrt(sum(predict(fit206)-x)^2)/(N-2))
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
  plot(x, y, main = "Dubinin-Radushkevich Isotherm Model", xlab = "Ce",
       ylab = "Qe")
  lines(x,predict(fit206),lty=1, col="black", lwd=1)
  rsqq <- lm(Qe~predict(fit206))
  print(summary(rsqq))
}
