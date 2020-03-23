#' @title Baudu Isotherm Analysis Non-Linear Form
#' @description Baudu is a reduced form of the Langmuir Isotherm upon observation that the estimation of Langmuir coefficients, b and q, by measurement of tangents at different equilibrium concentrations shows that they are not constants in a broad range.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numericl value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"
#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Baudu isotherm
#' @examples bauduanalysis(c(.014, .063, .094, .2, .385, 1.15, 1.64), c(.072, .213, .282, .32, .325, .338, .344))
#' @author Roan Maeve E.Arilla
#' @author Daren Mae B. Imperial
#' @author C.C. Deocaris
#' @references Ayawei, N., Ebelegi, A. and Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, pp.1-12. doi: 10.1155/2017/3039817
#' @references Subramania, D. and Ramadoss, R. (2018). Adsorption of Chromium Using Blue Green Algae-Modeling and Application of Various Isotherms. International Journal of Chemical Technology, 10(1), pp.1-22.
#' @export
bauduanalysis <- function(Ce,Qe){
  x <- Ce
  y <- Qe
  data <- data.frame(x, y)
  print("NLS2 Analysis for Baudu Isotherm")
  fit200 <- Qe ~ ((Qmax*Ce*b)^(1+x1+y1))/((1+b*Ce)^(1+x1))
  start <- data.frame(Qmax= c(0,100), b= c(0,100), x1= c(0, 10), y1= c(-10,10))
  set.seed(511)
  fit201 <- nls2(fit200,
                data=data,
                start=start,
                control= nls.control(maxiter=1000, warnOnly=TRUE),
                algorithm= "plinear-random")
  print(summary(fit201))
  N <- nrow(na.omit(data))
  error <- function(y){
    pv  <- (predict(fit201))
    rmse<- (rmse(y,predict(fit201)))
    mae <- (mae(y,predict(fit201)))
    mse <- (mse(y,predict(fit201)))
    rae <- (rae(y,predict(fit201)))
    PAIC <- AIC(fit201)
    PBIC <- BIC(fit201)
    SE <-(sqrt(sum(predict(fit201)-x)^2)/(N-2))
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
  plot(Ce,Qe, main ="Baudu Isotherm Plot", xlab= "Ce", ylab= "Qe")
  lines(x, predict(fit201), col= "black")
  rsqq <- lm(Qe~predict(fit201))
  print(summary(rsqq))
}
