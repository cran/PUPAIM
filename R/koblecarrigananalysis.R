#' @title Koble-Carrigan Isotherm
#' @description It is three-parameter isotherm model equation that incorporates both Freundlich and Langmuir isotherms for representing equilibrium adsorption data. Koble-Corrigan isotherm model appeared to have advantages over both the Langmuir and Freundlich equations in that it expresses adsorption data over very wide ranges of pressures and temperatures.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Koble-Carrigan isotherm analysis
#' @examples koblecarrigananalysis(moringa$Ce, moringa$Qe)
#' @author Reinald L. Claudio
#' @author Verna Chin DR. Caparanga
#' @author C.C. Deocaris
#' @references Ayawei, N., Ebelegi, A. N., & Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 1-11. doi: 10.1155/2017/3039817
#' @references Koble, R., & Corrigan, T., "Adsorption isotherms for pure hydrocarbons," Industrial and Engineering Chemistry, vol. 44, no. 2, pp. 383-387, 1952.
#' @export
koblecarrigananalysis <- function(Ce,Qe){
  x <- Ce
  y <- Qe
  data <- data.frame(Ce, Qe)
  print("NLS2 Analysis for Koble-Carrigan Isotherm")
  mod245 <- Qe ~ ((Ak*Bk*Ce^p)/(1 + Bk*Ce))
  start <- data.frame(Ak = c(0, 1000), Bk = c(0, 1000), p = c(0, 1))
  set.seed(511)
  fit246 <- nls2(mod245, start = start, control = nls.control(maxiter =50 , warnOnly = TRUE),
            algorithm = c("port"))
  print(summary(fit246))
  N <- nrow(na.omit(data))
  error <- function(y){
    pv  <- (predict(fit246))
    rmse<- (rmse(y,predict(fit246)))
    mae <- (mae(y,predict(fit246)))
    mse <- (mse(y,predict(fit246)))
    rae <- (rae(y,predict(fit246)))
    PAIC <- AIC(fit246)
    PBIC <- BIC(fit246)
    SE <-(sqrt(sum(predict(fit246)-x)^2)/(N-2))
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
  plot(x, y, main ="Koble-Carrigan Isotherm Plot", xlab= "Ce", ylab= "Qe")
  lines(smooth.spline(x, predict(fit246)), col = "black")
  rsqq <- lm(Qe~predict(fit246))
  print(summary(rsqq))
}
