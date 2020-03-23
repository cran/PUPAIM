#' @title Fritz-Schluender Isotherm Analysis Non-Linear Form
#' @description An empirical equation which can fit a wide range of experimental results because of the large number of coefficients in the isotherm.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "plot"
#' @importFrom graphics "abline"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Fritz-Schlunder isotherm
#' @examples fritzanalysis(moringa$Ce, moringa$Qe)
#' @author Diego Calamba
#' @author Rain Grace R. Reyes
#' @author C.C. Deocaris
#' @references Aljeboree, A. M., Alshirifi, A. N., & Alkaim, A. F. (2017). Kinetics and equilibrium study for the adsorption of textile dyes on coconut shell activated carbon. Arabian Journal of Chemistry, 10, S3381-S3393. doi: 10.1016/j.arabjc.2014.01.020
#' @export
fritzanalysis <- function (Ce, Qe){
  x <- Ce
  y <- Qe
  data <- data.frame(Ce, Qe)
  print("NLS2 Analysis for Fritz-Schlunder Analysis")
  mod221 <- Qe ~ (C*Ce^a)/(1+D*Ce^b)
  start <- data.frame(C = c(0, 1000),
                      D = c(0, 1000),
                      a = c(0, 1),
                      b = c(0, 1))
  set.seed(511)
  fit222 <- nls2(mod221, start = start, control = nls.control(maxiter =1000 , warnOnly = TRUE),
               algorithm = "plinear-random")
  pars <- as.vector(coefficients(fit222))
  pars_C <- pars[1L]; pars_D <- pars[2L]; pars_a <- pars[3L]; pars_b <- pars[4L]; pars_lin <- pars[5L]
  start1 <- data.frame(C = pars_C*pars_lin,
                       D = pars_D,
                       a = pars_a,
                       b = pars_b)
  fit223 <- nls2(mod221,
               start = start1,
               algorithm = "brute-force",
               control=list(maxiter=100))
  print(summary(fit223))
  N <- nrow(na.omit(data))
  error <- function(y){
    pv  <- (predict(fit223))
    rmse<- (rmse(y,predict(fit223)))
    mae <- (mae(y,predict(fit223)))
    mse <- (mse(y,predict(fit223)))
    rae <- (rae(y,predict(fit223)))
    PAIC <- AIC(fit223)
    PBIC <- BIC(fit223)
    SE <-(sqrt(sum(predict(fit223)-x)^2)/(N-2))
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
    plot(x, y, main ="Fritz-Schlunder Four Parameter Isotherm Non-Linear Plot", xlab= "Ce", ylab= "Qe")
  lines(smooth.spline(x,predict(fit223)), col = "black")
  rsqq <- lm(Qe~predict(fit223))
  print(summary(rsqq))
}
