#' @title Flory-Huggins Isotherm Analysis Non-Linear Form
#' @description Flory-Huggins isotherm is a two-parameter isotherm that describes the degree of surface coverage characteristics of the adsorbate on the adsorbent. This isotherm model can express the feasibility and spontaneity of an adsorption process.
#' @param Ce is equal to Co which is the numeric value for the iniztial concentration
#' @param Qe is equal to theta which is the degree of surface coverage
#' @importFrom graphics "plot"
#' @importFrom graphics "abline"
#' @importFrom nls2 "nls2"
#' @importFrom stats "predict" "lm" "AIC" "BIC"
#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Flory-Huggins isotherm
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples fhanalysis(Ce, Qe)
#' @author Carmela L. Barbacena
#' @author C.C. Deocaris
#' @references Ayawei, N., et al., (5 September 2017). Modeling and Interpretation of Adsorption Isotherms. Retrieved from: https:// www.hindawi.com/journals/jchem/2017/3039817/
#' @references Tsamo, C., et al., (27 November 2019). One-, Two-, and Three- Parameter Isotherms, Kinetics, and Thermodynamic Evaluation of Co
#' @export
fhanalysis <- function(Ce, Qe){
  x <- 1-Qe
  y <- Qe/Ce
  dat <- data.frame(x,y)
  nx<- nrow(na.omit(dat))
  mod213 <- (y ~ kfh*((x)^n))
  start <- data.frame(kfh = c(0, 100), n = c(0, 100))
  set.seed(511)
  fit214 <- nls2(mod213, start = start, control = nls.control(maxiter = 50, warnOnly = TRUE), algorithm = "brute-force")
  print("Flory Huggins Parameters")
  print(summary(fit214))
  pv <- predict(fit214)
  errors <- function(y) {
    rmse <-rmse(y, predict(fit214))
    mae <- mae(y, predict(fit214))
    mse <- mse(y, predict(fit214))
    rae <- rae(y, predict(fit214))
    PAIC <- AIC(fit214)
    PBIC <- BIC(fit214)
    SE <-(sqrt(sum(predict(fit214)-x)^2)/(nx-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Predicted Values" = pv,"relative mean squared error" = rmse, "mean absolute error" = mae, "mean squared error" = mse, "relative absolute error" = rae,
         "AIC" = PAIC,  "BIC" = PBIC, "Standard Error Estimate" = SE)
  }
  f <- errors(y)
  print(f)
  plot(x, y, main = "Flory-Huggins Isotherm Model", xlab = "1-Qe",
       ylab = "Qe/Ce")
  lines(x,predict(fit214),lty=1,col="black",lwd=1)
  rsqq <- lm(Qe~predict(fit214))
  print(summary(rsqq))
}
