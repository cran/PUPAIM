#' @title Toth Isotherm Analysis Non-Linear Form
#' @description Another empirical modification of the Langmuir equation with the aim of reducing the error between experimental data and predicted value of equilibrium data.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the fractional coverage
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return The non linear regression and the parameters for the Toth isotherm analysis
#' @examples tothanalysis(c(1,2,3,4,5),c(1,2,3,4,5))
#' @author Kim Zyrell P. Zagala
#' @author C.C. Deocaris
#' @references Gutierrez, L.G., et.al(2018, September 20), Kinetic and Equilibrium Study of the Absorption of CO2 in Ultramicropores of Resorcinol-Formaldehyde Aerogels obtained in Acidic and Basic Medium. Retrieved from: doi:10.3390/c4040052
#' @references Ayawei, N. (2017, September 05). Medelling an Interpretation of Adsorption Isotherm. Retrieved from: https://www.hidawi.com/journals/jchem/2017/3039817
#' @export
tothanalysis <- function(Ce, Qe){
  x <- Ce
  y <- Qe
  fit265 <- (y ~ ((Qt)*(Kt)*(x))/((1+((Kt)*(x))^(Mt))^(1/(Mt))))
  start <- data.frame(Qt = c(0, 100), Kt = c(0, 100), Mt = c(0, 100))
  set.seed(511)
  suppressWarnings(fit266 <- nls2(fit265, start = start, control = nls.control(maxiter = 100, warnOnly = TRUE), algorithm = "plinear-random"))
  print(summary(fit266))
  N <- nrow(na.omit(data))
  error <- function(y){
    pv  <- (predict(fit266))
    rmse<- (rmse(y,predict(fit266)))
    mae <- (mae(y,predict(fit266)))
    mse <- (mse(y,predict(fit266)))
    rae <- (rae(y,predict(fit266)))
    PAIC <- AIC(fit266)
    PBIC <- BIC(fit266)
    SE <-(sqrt(sum(predict(fit266)-x)^2)/(N-2))
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
    plot(x, y, main = "Toth Isotherm Analysis", xlab = "Ce", ylab = "Qe")
    lines(x, predict(fit266), col = "black")
    rsqq <- lm(Qe~predict(fit266))
    print(summary(rsqq))
}
