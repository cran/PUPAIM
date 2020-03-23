#' @title Hill Isotherm Analysis Non-Linear Form
#' @description Hill isotherm model shows the connection of different species of homogeneous surfaces. It assumes that the adsorption is a cooperative phenomenon, with a ligand binding activity at one  part of a macromolecule that may affect the different binding sites of that same macromolecule.
#' @param Ce equilibrium capacity
#' @param Qe adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @examples Ce <- c(0.001, 0.0026, 0.0125, 0.031, 0.056)
#' @examples Qe <- c(0.02, 0.072, 0.146, 0.15, 0.151)
#' @examples hillanalysis(Ce, Qe)
#' @author Amiela D. Suerte
#' @author Carl Luis P. Flestado
#' @author C.C. Deocaris
#' @references Tanzifi, M., Karimipour, K., Hoseini, S., Ali, I.(2017). Artificial neural network optimization for methyl orange adsorption onto polyaniline nano-adsorbent: Kinetic, isotherms and thermodynamics. Journal of Molecular Liquids, 2017, p.11. DOI:10.1016/j.molliq.2017.08.122
#' @references Saadi, R., Saadi, z., Fazaeli, R., Fard, N.E.(2015). Monolayer and multilayer adsorption models for sorption aqueous media. Korean Journal of Chemical Engineering, 2015, p.5. DOI: 10.007/s11814-015-0053-7
#' @references Larimi, S.G., Moghadamnia, A.A., Najafpour, G.(2016). Kinetics and isotherm studies of the Immobilized Lipase on Chitosan Support. International Journal of Engineering, 2026, p.12. DOI: 10.5829/idosi.ije.2016.29.10a.01
#' @export
hillanalysis <- function(Ce, Qe){
  x <- Ce
  y <- Qe
  dat <- data.frame(x,y)
  n<- nrow(na.omit(dat))
  mod231 <- (y ~ ((qh*(x^nh))/(kd+x^(nh))))
  print("Hill Analysis")
  start <- data.frame(qh = c(0, 100), nh = c(0, 2), kd = c(0, 10))
  set.seed(511)
  fit232 <- nls2(mod231, start = start, control = nls.control(maxiter = 100, warnOnly = TRUE), algorithm = "plinear-random")
  print(summary(fit232))
  error <- function(y){
    pv  <- (predict(fit232))
    rmse<- (rmse(y,predict(fit232)))
    mae <- (mae(y,predict(fit232)))
    mse <- (mse(y,predict(fit232)))
    rae <- (rae(y,predict(fit232)))
    PAIC <- AIC(fit232)
    PBIC <- BIC(fit232)
    SE <-(sqrt(sum(predict(fit232)-x)^2)/(n-2))
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
  plot(x, y, main = "Hill Isotherm", xlab = "Ce", ylab = "Qe")
  lines(smooth.spline(x, predict(fit232)))
  rsqq <- lm(Qe~predict(fit232))
  print(summary(rsqq))
}
