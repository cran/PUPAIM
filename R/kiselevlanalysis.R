#' @title Kiselev Isotherm Non linear Analysis
#' @description It is also known as localized monomolecular layer model and is only valid for surface coverage theta > 0.68.
#' @param theta the numerical value for surface coverage
#' @param Ce the numerical value for equilibrium capacity
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mse" "mae" "rae"
#' @return the linear regression and the parameters for the Kiselev isotherm analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples kiselevanalysis(Ce, Qe)
#' @author Ashley Quebrado
#' @author C. C. Deocaris
#' @export
kiselevanalysis <- function(theta, Ce){
    x <- theta
    y <- Ce
    data <- data.frame(x, y)
    print("NLS2 Analysis for Kiselev Isotherm Analysis")
    mod242 <- (Ce)~ (theta/(K1*(1-theta) * (1 + Ki*theta)))
    start <- data.frame(K1 = c(-1000,1000), Ki = c(-1000,1000))
    set.seed(511)
    suppressWarnings(fit243 <- nls2(mod242, start = start, control = nls.control(maxiter= 100, warnOnly = TRUE), algorithm = "port"))
    print(summary(fit243))
    N <- nrow(na.omit(data))
    error <- function(Ce){
      pv  <- (predict(fit243))
      rmse<- (rmse(Ce,predict(fit243)))
      mae <- (mae(Ce,predict(fit243)))
      mse <- (mse(Ce,predict(fit243)))
      rae <- (rae(Ce,predict(fit243)))
      PAIC <- AIC(fit243)
      PBIC <- BIC(fit243)
      SE <-(sqrt(sum(predict(fit243)-theta)^2)/(N-2))
      colnames(Ce) <- rownames(Ce) <- colnames(Ce)
      list("Predicted Values"           = pv,
           "Relative Mean Square Error" = rmse,
           "Mean Absolute Error"        = mae,
           "Mean Squared Error"         = mse,
           "Relative Absolute Error"    = rae,
           "AIC"                        = PAIC,
           "BIC"                        = PBIC,
           "Standard Error Estimate"    = SE)
    }
    e <- error(Ce)
    print(e)
    plot(x, y, main = "Kiselev Isotherm", xlab = "theta",
         ylab = "Ce")
    lines(x, predict(fit243), col = "black")
    rsqq <- lm(Ce~predict(fit243))
    print(summary(rsqq))
}
