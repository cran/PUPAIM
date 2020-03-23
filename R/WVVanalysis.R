#' @title Weber Van-Vliet Isotherm Analysis Non-linear Form
#' @description it provides an excellent description of data patterns for a broad range of systems. This model is suitable for batch rate and fixed-bed modelling procedures as it gives a direct parameter evaluation.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for Weber-Van-Vliet  Isotherm Analysis
#' @examples webervvanalysis(moringa$Qe,moringa$Ce)
#' @author Jeann M. Bumatay
#' @author Leslie F. Bautista
#' @author Chester C. Deocaris
#' @references Van Vliet, B.M., Weber Jr., Hozumi, H.. (1979). Modeling and prediction of specific compound adsorption by activated carbon and synthetic adsorbents. Water Research Vol.14, pp. 1719 to 1728. https://doi.org/10.1016/0043-1354(80)90107-4
#' @references
#' @export
webervvanalysis<- function(Qe,Ce) {
  x <- Qe
  y <- Ce
  data <- data.frame(Qe, Ce)
  n<- nrow(na.omit(data))
  print("NLS2 Analysis for Weber Van-Vliet Isotherm")
  fit267 <- (Ce ~ (P * Qe^(R*(Qe^s)+t)))
  start <- data.frame(P = c(1e-6, 10),
                      R = c(1e-7, 0.999),
                      s = c(1e-2, 3),
                      t = c(0.1, 5))
  set.seed(511)
  suppressWarnings(fit268 <- nls2(fit267,
                                 start = start,
                                 control = nls.control(maxiter = 50, warnOnly = TRUE),
                                 algorithm = c("plinear-random")))
  print(summary(fit268))
  error <- function(y){
    pv  <- (predict(fit268))
    rmse<- (rmse(y,predict(fit268)))
    mae <- (mae(y,predict(fit268)))
    mse <- (mse(y,predict(fit268)))
    rae <- (rae(y,predict(fit268)))
    PAIC <- AIC(fit268)
    PBIC <- BIC(fit268)
    SE <-(sqrt(sum(predict(fit268)-x)^2)/(n-2))
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
  plot(x, y, main ="Weber Van-Vliet Isotherm Plot", xlab= "Qe", ylab= "Ce")
  lines(x, predict(fit268), col = "black")
  rsqq <- lm(Ce~predict(fit268))
  print(summary(rsqq))
}
