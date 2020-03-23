#' @title Radke-Prausniiz Isotherm Analysis Nonlinear Form
#' @description Nonlinear form of the Radke-Prausniiz Equation, the original form. The Radke-Prausnitz isotherm model has several important properties which makes it more preferred in most adsorption systems at low adsorbate concentration.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Raudke-Prausniiz isotherm
#' @examples raudkepanalysis(moringa$Ce,moringa$Qe)
#' @author Princess Joyce DL Reyes
#' @author Neil Ross S. Alayon
#' @author C.C. Deocaris
#' @references Newton, A., & Donbebe. (2017, September 5). Modelling and Interpretation of Adsorption Isotherms. Retrieved from https://doi.org/10.1155/2017/3039817
#' @references Khalid, A., Kazmi, M. et.al, (2015). Kinetic & Equilibrium Modelling of Copper Biosorption. Retrieved from http://journals.pu.edu.pk/journals/index.php/jfet/article/view/527
#' @export
raudkepanalysis <- function(Ce, Qe){
  x <- Ce
  y <- Qe
  dat <- data.frame(x,y)
  n<- nrow(na.omit(dat))
  mod254 <- (y ~ (Krp*Frp)/((Krp*(x^(1-Nrp)))+Frp))
  print("Raudke Analysis")
  start <- data.frame(Krp = c(-100, 100), Frp = c(0, 1000), Nrp = c(0, 1)) ##sure ka?
  set.seed(511)
  fit255 <- nls2(mod254, start = start, control = nls.control(maxiter = 100 , warnOnly = TRUE), algorithm = "port")
print(summary(fit255))
error <- function(y){
  pv  <- (predict(fit255))
  rmse<- (rmse(y,predict(fit255)))
  mae <- (mae(y,predict(fit255)))
  mse <- (mse(y,predict(fit255)))
  rae <- (rae(y,predict(fit255)))
  PAIC <- AIC(fit255)
  PBIC <- BIC(fit255)
  SE <-(sqrt(sum(predict(fit255)-x)^2)/(n-2))
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
  plot(x, y, main = "Raudke-Prausniiz Isotherm", xlab = "Ce",
       ylab = "Qe")
   lines(smooth.spline(Ce, predict(fit255)), col = "black")
   rsqq <- lm(Qe~predict(fit255))
   print(summary(rsqq))
}
