#' @title Freundlich Isotherm Analysis Linear Form
#' @description This isotherm model predicts the logarithmic relationship of equilibrium capacity and  adsorbed capacity (Gessner and Hasan, 1987). Furthermore, this model gives an equation which defines the surface heterogeneity and the exponential distribution of active sites (Ayawei, et al., 2017).
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "plot" "abline"

#' @importFrom Metrics "rmse" "mse" "mae" "rae"
#' @return the linear regression and the parameters for the Freundlich isotherm
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples freundlich.LM (Ce, Qe)
#' @author Carl Luis P. Flestado
#' @author Chester C. Deocaris
#' @references Ayawei, N., Ebelegi, A. N., & Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 1–11. doi: 10.1155/2017/3039817
#' @references Gessner, P. K., & Hasan, M. M. (1987). Freundlich and Langmuir Isotherms as Models for the Adsorption of Toxicants on Activated Charcoal. Journal of Pharmaceutical Sciences, 76(4), 319–327. doi: 10.1002/jps.2600760412
#' @export
freundlich.LM <- function (Ce, Qe)
{
  x <- log10(Ce)
  y <- log10(Qe)
  dat <- data.frame(x,y)
  N<- nrow(na.omit(dat))
  rhs <- function(x, Kf, n) {
    log10(Kf)+(1/n)*log10(x)
  }
  fit220 <- lm(y~x)
  print("Freundlich Isotherm Analysis")
  print(summary(fit220))
  a <- (summary(fit220))
  b <- a$coefficients[2]
  c <- a$coefficients[1]

  param.ads <- function(y){
    Kf <- 10^(c)
    n <- 1/b
    colnames(y) <- rownames(y) <- colnames(y)
    print("Freundlich Parameters")
    list("Kf" = Kf,
         "n" = n)
  }
  w <- param.ads(y)
  print(w)
  error <- function(y){
    pv  <- (predict(fit220))
    rmse<- (rmse(y,predict(fit220)))
    mae <- (mae(y,predict(fit220)))
    mse <- (mse(y,predict(fit220)))
    rae <- (rae(y,predict(fit220)))
    PAIC <- AIC(fit220)
    PBIC <- BIC(fit220)
    SE <-(sqrt(sum(predict(fit220)-x)^2)/(N-2))
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
  plot(x, y, main = "Freundlich Linear Isotherm Model", xlab = "log10(Ce)",
       ylab = "log10(Qe)")
  abline(fit220, col="black")
}
