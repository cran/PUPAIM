#' @title Sips Isotherm Analysis Non linear Form
#' @description This model is suitable for predicting adsorption on heterogeneous surfaces, thereby avoiding the limitation of increased adsorbate concentration normally associated with the Freundlich model.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Sips isotherm
#' @examples sipsanalysis(moringa$Ce,moringa$Qe)
#' @references N'diaye, A., Bollahi, M.,Kankou, M. (2019). Sorption of paracetamol from aqueous solution using groundnut shell as a low cost sorbent. J. Mater. Environ. Sci., 2019, Vol.10, Issue 6, 553-562.
#' @references Nethaji, S.,Sivasamy, A., Mandal, A. B. (2012). Adsorption isotherms, kinetics and mechanism for the adsorption of cationic and anionic dyes onto carbonaceous particles prepared from Juglans regia shell biomass. Int. J. Environ. Sci. Technol. (2013)10:231-242. doi: 10.1007/s13762-012-0112-0
#' @export
sipsanalysis <- function(Ce, Qe){
  x <- Ce
  y <- Qe
  dat <- data.frame(x,y)
  n<- nrow(na.omit(dat))
  fit259 <- (y ~ ((Qm*Ks*x^n)/(1 + Ks*x^n)))
  print("Sips Analysis")
  start <- data.frame(Qm = c(0, 100), Ks = c(0, 100), n = c(0, 10))
  set.seed(511)
  suppressWarnings(fit260 <- nls2(fit259, start = start, control = nls.control(maxiter = 30, warnOnly = TRUE), algorithm = "port"))
  print(summary(fit260))
  predict(fit260)
  error <- function(y){
    pv  <- (predict(fit260))
    rmse<- (rmse(y,predict(fit260)))
    mae <- (mae(y,predict(fit260)))
    mse <- (mse(y,predict(fit260)))
    rae <- (rae(y,predict(fit260)))
    PAIC <- AIC(fit260)
    PBIC <- BIC(fit260)
    SE <-(sqrt(sum(predict(fit260)-x)^2)/(n-2))
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
  nplot <- function(Ce, Qe) {
    x <- (Ce); y <- (Qe)
    plot(x, y, main = "Sips Isotherm", xlab = "Ce", ylab = "Qe")
    lines(smooth.spline(Ce, predict(fit260)), col = "black")
  }
  nplot(Ce, Qe)
  rsqq <- lm(Qe~predict(fit260))
  print(summary(rsqq))
}
