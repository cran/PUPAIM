#' @title Redlich-Peterson Non-Linear Form
#' @description It is used for three parameter adsorption and its a combination of Langmuir and Freundlich Isotherms
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the linear regression and the parameters for the Redlich-Peterson isotherm
#' @examples redlichpanalysis(moringa$Ce,moringa$Qe)
#' @author John Carlo F. Panganiban
#' @author C.C. Deocaris
#' @references Ayawei, N., Ebelegi, A. N., & Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 1.11. doi: 10.1155/2017/3039817
#' @references Wu, F.-C., Liu, B.-L., Wu, K.-T., & Tseng, R.-L. (2010). A new linear form analysis of Redlich?Peterson isotherm equation for the adsorptions of dyes. Chemical Engineering Journal, 162(1), 21.27. doi: 10.1016/j.cej.2010.03.006
#' @export
redlichpanalysis <- function(Ce, Qe){
    x <- Ce
    y <- Qe
    data <- data.frame(x, y)
    n<- nrow(na.omit(data))
    print("NLS2 Analysis for Redlich Peterson Isotherm Model")
    fit257 <- (Qe ~ (Arp*Ce)/(1+(Krp*(Ce^b))))
    start <- data.frame(Arp = c(1, 100), Krp = c(1, 100), b = c(0, 1))
    set.seed(511)
    fit258 <- nls2(fit257, start = start, control = nls.control(maxiter = 50, warnOnly = TRUE), algorithm = "port")
    print(summary(fit258))
    error <- function(y){
      pv  <- (predict(fit258))
      rmse<- (rmse(y,predict(fit258)))
      mae <- (mae(y,predict(fit258)))
      mse <- (mse(y,predict(fit258)))
      rae <- (rae(y,predict(fit258)))
      PAIC <- AIC(fit258)
      PBIC <- BIC(fit258)
      SE <-(sqrt(sum(predict(fit258)-x)^2)/(n-2))
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
    plot(x, y, main = "Redlich Peterson Isotherm", xlab = "Ce", ylab = "Qe")
    lines(x, predict(fit258), col = "black")
    rsqq <- lm(predict(fit258)~Qe)
    print(summary(rsqq))
  }
