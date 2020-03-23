#' @title Freundlich Isotherm Analysis Non-Linear Form
#' @description This isotherm model is an empirical model applicable to diluted solutions adsorption processes (Gessner and Hasan, 1987). Furthermore, this model gives an equation which defines the surface heterogeneity and the exponential distribution of active sites (Ayawei, et al., 2017).
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Freundlich isotherm
#' @examples  Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples  Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @author Carl Luis P. Flestado
#' @author C.C. Deocaris
#' @references Ayawei, N., Ebelegi, A. N., & Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 1–11. doi: 10.1155/2017/3039817
#' @references Gessner, P. K., & Hasan, M. M. (1987). Freundlich and Langmuir Isotherms as Models for the Adsorption of Toxicants on Activated Charcoal. Journal of Pharmaceutical Sciences, 76(4), 319–327. doi: 10.1002/jps.2600760412
#' @export
freundlichanalysis <- function(Ce,Qe){
  x <- Ce
  y <- Qe
  fxn <- y ~ Kf*x^(1/n)
  dat <- data.frame(x,y)
  n<- nrow(na.omit(dat))
  grd <- data.frame(Kf=c(0,1000),
                    n=c(0,10))
  set.seed(511)
  fit1 <- nls2(fxn,
               start = grd,
               algorithm = "plinear-random",
               control=list(maxiter=1000))
  pars <- as.vector(coefficients(fit1))
  pars_Kf <- pars[1L]; pars_n <- pars[2L]; pars_lin <- pars[3L]
  grd1 <- data.frame(Kf=pars_Kf*pars_lin,
                     n=pars_n)
  fit2 <- nls2(fxn,
               start = grd1,
               algorithm = "brute-force",
               control=list(maxiter=100))
  print(summary(fit2))
  error <- function(y){
    pv  <- (predict(fit2))
    rmse<- (rmse(y,predict(fit2)))
    mae <- (mae(y,predict(fit2)))
    mse <- (mse(y,predict(fit2)))
    rae <- (rae(y,predict(fit2)))
    PAIC <- AIC(fit2)
    PBIC <- BIC(fit2)
    SE <-(sqrt(sum(predict(fit2)-x)^2)/(n-2))
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
  plot(dat$x,dat$y,
       main = "Freundlich Isotherm Non-Linear Plot", xlab = "Ce", ylab= "Qe")
  lines(smooth.spline(x, predict(fit2)), col = "black")
  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))
}
