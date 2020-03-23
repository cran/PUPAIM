#' @title Five Parameter Isotherm Analysis
#' @description A five-parameter empirical model that is capable of simulating the model variations that applies over a wide range of equilibrium data. Its increased parameters provides a more accurate non linear regression, better than two, three or four parameters. (Subramanyam and Ashutosh, 2009)
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"
#' @importFrom stats "predict" "lm" "AIC" "BIC" "nls.control" "na.omit" "smooth.spline" "coefficients"
#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the non linear regression and the parameters for the five parameter non-linear isotherm analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples fiveparamanalysis(Ce, Qe)
#' @author Jeff Ryan S. Magalong
#' @author Joshua Dela Cruz
#' @author Chester C. Deocaris
#' @references Fritz, W., & Schluender, E.U. (1974). Simultaneous adsorption equilibria of organic solutes in dilute aqueous solutions on activated carbon. Chemical Engineering Science, 29(5), 1279-1282. <doi: 10.1016/0009-2509(74)80128-4>
#' @references @references Subramanyam, B., & Das, A. (2009). Study of the adsorption of phenol by two soils based on kinetic and isotherm modeling analyses. Desalination, 249(3), 914-921. <doi:10.1016/j.desal.2009.05.020>
#' @export
fiveparamanalysis <- function(Ce, Qe){
  x <- Ce
  y <- Qe
  mod210 <- (y ~ (Qmax* K1*x^alpha)/(1+ K2*x^beta))
  dat <- data.frame(x,y)
  n<- nrow(na.omit(dat))
  grd <- data.frame(Qmax=c(0,1000),
                    K1=c(0,10),
                    K2=c(0,100),
                    alpha= c(0,1),
                    beta= c(0,1))
  set.seed(511)
  fit211 <- nls2(mod210,
               start = grd,
               algorithm = "plinear-random",
               control=list(maxiter=1000))
  pars <- as.vector(coefficients(fit211))
  pars_Qmax <- pars[1L]; pars_K1 <- pars[2L]; pars_K2 <- pars[3L]; pars_alpha <- pars[4L]; pars_beta <- pars[5L]; pars_lin <- pars[6L]
  grd1 <- data.frame(Qmax=pars_Qmax*pars_lin,
                     K1=pars_K1,
                     K2=pars_K2,
                     alpha=pars_alpha,
                     beta= pars_beta)
  fit212 <- nls2(mod210,
               start = grd1,
               algorithm = "brute-force",
               control=list(maxiter=1000))
  print(summary(fit212))
  error <- function(y){
    pv  <- (predict(fit212))
    rmse<- (rmse(y,predict(fit212)))
    mae <- (mae(y,predict(fit212)))
    mse <- (mse(y,predict(fit212)))
    rae <- (rae(y,predict(fit212)))
    PAIC <- AIC(fit212)
    PBIC <- BIC(fit212)
    SE <-(sqrt(sum(predict(fit212)-x)^2)/(n-2))
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
       main = "Fritz-Schluender Five Parameter Isotherm Non-Linear Plot", xlab = "Ce", ylab= "Qe")
  lines(smooth.spline(x, predict(fit212)), col = "black")
  rsqq <- lm(Qe~predict(fit212))
  print(summary(rsqq))
}
