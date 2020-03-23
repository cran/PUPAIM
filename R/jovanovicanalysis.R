#' @title Jovanovic Isotherm Analysis Non-Linear Form
#' @description It is predicated on the assumptions contained in the Langmuir model, but in addition the possibility of some mechanical contacts between the adsorbate and adsorbent
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorpted capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#'@return the nonlinear regression and the parameters for the Jovanovic isotherm
#'@examples jovanovicanalysis(moringa$Ce, moringa$Qe)
#' @author Christian Josuah F. Maylas
#' @author C.C. Deocaris
#' @references:  Saadi, R., Saadi, Z., Fazaeli, R., Fard, N. E. (2015). Monolayer and multilayer adsorption isotherm models for sorption from aqueous media. Korean J. Chem. Eng., 32(5), 787-799 (2015) DOI: 10.1007/s11814-015-0053-7
#' @references:  Vargas, A., Cazetta, A., Kunita, M., Silva, T., Almeida V. (2011). Adsorption of methylene blue on activated carbon produced from ?amboyant pods (Delonix regia): Study of adsorption isotherms and kinetic models. Chemical Engineering Journal 168 (2011) 722-730
#' @export
jovanovicanalysis <- function(Ce,Qe){
  x <- Ce
  y <- Qe
  fxn <- y ~ Qm*(1-exp(-kj*x))
  dat <- data.frame(x,y)
  n<- nrow(na.omit(dat))
  grd <- data.frame(Qm=c(0,1000),
                    kj=c(0,10))
  fit1 <- nls2(fxn,
               start = grd,
               algorithm = "plinear-random",
               control=list(maxiter=1000))
  pars <- as.vector(coefficients(fit1))
  pars_Qm <- pars[1L]; pars_kj <- pars[2L]; pars_lin <- pars[3L]
    grd1 <- data.frame(Qm=pars_Qm*pars_lin,
                     kj=pars_kj)
  fit2 <- nls2(fxn,
               start = grd1,
               algorithm = "brute-force",
               control=list(maxiter=100))
  print(summary(fit2))
  error <- function(y){
    pv  <- predict(fit2)
    rmse<- rmse(y,predict(fit2))
    mae <- mae(y,predict(fit2))
    mse <- mse(y,predict(fit2))
    rae <- rae(y,predict(fit2))
    PAIC <- AIC(fit2)
    PBIC <- BIC(fit2)
    SE <-(sqrt(sum(predict(fit2)-x)^2)/(n-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Predicted Values"           = pv,
         "Relative Mean Square Error" = rmse,
         "Mean Absolute Error"        = mae,
         "Mean Squared Error"         = mse,
         "Relative Absolute Error"    = rae,
         "AIC"                        = PAIC ,
         "BIC"                        = PBIC,
         "Standard Error Estimate"    = SE)
  }
  e <- error(y)
  print(e)
  plot(dat$x,dat$y,
       main = "Jovanovic Isotherm Non-Linear Plot", xlab = "Ce", ylab= "Qe")
  lines(smooth.spline(x, predict(fit2)), col = "black")
  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))
}
