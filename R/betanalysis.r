#' @title BET Isotherm Analysis Non-Linear Form
#' @description An isotherm that takes account of the possibility that the monolayer in the Langmuir adsorption isotherm can act as a substrate for further adsorption.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "plot" "abline" "lines"
#' @importFrom nls2 "nls2"
#' @importFrom stats "predict" "lm" "AIC" "BIC" "nls.control" "na.omit"
#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the non-linear regression, errors, and the parameter values for the BET Isotherm Model
#' @examples BETanalysis(c(.014, .063, .094, .2, .385, 1.15, 1.64), c(.072, .213, .282, .32, .325, .338, .344))
#' @author Geraldine N. Aliman
#' @author Cherry Rose E. Olivar
#' @author C.C. Deocaris
#' @references Shipu, S.(2014, April 30). Bet Isotherm. Retrived form https://www.slideshare.net/mobile/Sourav44Shipu/bet-isotherm
#' @export
BETanalysis <- function(Ce,Qe) {
  x <- Ce
  y <- Qe
  data <- data.frame(Qe, Ce)
  print("NLS2 Analysis for BET Isotherm Analysis")
  fit202 <- Qe ~ ((Qmax*Kb*Ce)/((Cs-Ce)*(1+(Kb-1)*(Ce/Cs))))
  start <- data.frame(Qmax = c(-1, 100),
                      Kb = c(0, 100),
                      Cs = c(-1, 1))
  set.seed(511)
  suppressWarnings(fit203 <- nls2(fit202,
                start = start,
                control = nls.control(maxiter =100 , warnOnly = TRUE),
                algorithm = "port"))
  print(summary(fit203))
  N <- nrow(na.omit(data))
  error <- function(y){
    pv  <- (predict(fit203))
    rmse<- (rmse(y,predict(fit203)))
    mae <- (mae(y,predict(fit203)))
    mse <- (mse(y,predict(fit203)))
    rae <- (rae(y,predict(fit203)))
    PAIC <- AIC(fit203)
    PBIC <- BIC(fit203)
    SE <-(sqrt(sum(predict(fit203)-x)^2)/(N-2))
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
  plot(x, y, main = "BET Analysis", xlab = "Ce", ylab = "Qe")
  lines(x, predict(fit203), col = "black")
  rsqq <- lm(Qe~predict(fit203))
  print(summary(rsqq))
}
