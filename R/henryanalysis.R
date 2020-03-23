#' @title Henry's Isotherm
#' @description It describes the appropriate fit to the adsorption of adsorbate at relatively low concentrations such that all adsorbate molecules are secluded from their nearest neighbours.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the linear regression and the parameter for the henry isotherm analysis
#' @importFrom graphics "plot" "abline"

#' @importFrom minpack.lm "nlsLM"
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples henryanalysis (Ce, Qe)
#' @author De Osio, Lloyd P.
#' @references Ayawei, N., Ebelegi, A. N., and Wankasi, D. (2017). Review Article: Modelling and interpretation of Adsorption isotherms. HINDAWI: Journal of Chemistry. doi: 10.1155/2017/3039817
#' @references Bayuo, J. Pelig-Ba, K.B., and Abukari, M.A. (2018). Isotherm modelling of lead (II) adsorption from aqueous solution using groundnut shell as a low-cost adsorbent. IOSR Journal of Appplied chemistry doi:10.9790/5736-1111011823
#' @export
henryanalysis <- function (Ce, Qe)
{
  x <- Ce
  y <- Qe
  dat <- data.frame(x,y)
  N<- nrow(na.omit(dat))
  rhs <- function(x, K) {
    (K*x)
  }
  fit230 <- lm(y ~ x)
  print("Henry's Isotherm Analysis")
  print(summary(fit230))
  a <- (summary(fit230))
  b <- a$coefficients[2]
  param.ads <- function(y){
    qmax <- 1/a$coefficients[1]
    b <- (a$coefficients[2]/qmax)
    colnames(y) <- rownames(y) <- colnames(y)
    print("Henry's Isotherm Parameter")
    list("KH" = b)
  }
  w <- param.ads(y)
  print(w)
  error <- function(y){
    pv  <- (predict(fit230))
    rmse<- (rmse(y,predict(fit230)))
    mae <- (mae(y,predict(fit230)))
    mse <- (mse(y,predict(fit230)))
    rae <- (rae(y,predict(fit230)))
    PAIC <- AIC(fit230)
    PBIC <- BIC(fit230)
    SE <-(sqrt(sum(predict(fit230)-x)^2)/(N-2))
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
  plot(x, y, main = "Henry's Isotherm Analysis", xlab = "Ce",
       ylab = "Qe")
  lines(x, predict(fit230))
}
