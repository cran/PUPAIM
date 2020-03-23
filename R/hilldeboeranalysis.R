#' @title Hill-Deboer Isotherm Non-linear Analysis
#' @description Hill-Deboer isotherm model describes as a case where there is mobile adsorption as well as lateral interaction among molecules.
#' @param theta the numerical value for the surface converage / fractional coverage
#' @param Ce the numerical value for the equilibrium capacity
#' @param t the numerical value for the temperature of the adsorption experimentation in Kelvin
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"
#' @importFrom stats "predict" "lm" "AIC" "BIC"
#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @examples Ce <- c(0.39, 0.74, 1.19, 1.63, 2.03, 2.51, 2.96, 3.46, 3.9, 4.35)
#' @examples theta <- c(0.051, 0.12, 0.14, 0.17, 0.21, 0.22, 0.24, 0.24, 0.27, 0.29)
#' @examples hilldeboeranalysis(theta, Ce, 298)
#' @author Benz L. Rivera
#' @author Jeff Ryan S. Magalong
#' @author C.C. Deocaris
#' @references Ayawei, N., Ebelegi, A.N., & Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 1-11. doi 10.115520173039817
#' @references Qing Shao and Carol K Hall (2016). Protein adsorption on nanoparticles model development using computer simulation. Journal of Physics Condensed Matter
#' @references Oualid Hamdaouia, Emmanuel Naffrechoux (2007). Modeling of adsorption isotherms of phenol andchlorophenols onto granular activated carbonPart I. Two-parameter models and equations allowing determination of thermodynamic parameters
#' @export
hilldeboeranalysis <- function(theta, Ce, t){
  x <- theta
  y <- Ce
  mod <- (y ~((x/(1-x))*exp((x/1-x)-(x*K2)/8.314*t))/K1)
  dat <- data.frame(x,y)
  n <- nrow(na.omit(dat))
  grd <- data.frame(K1=c(0,100),
                    K2=c(-10,100))
  set.seed(429)
  fit <- nls2(mod,
              start = grd,
              algorithm = "plinear-random",
              control = list(maxiter=1000),
              lower=c(0,-10),upper=c(100,100))
  pars <- as.vector(coefficients(fit))
  pars_K1<- pars[1L]; pars_K2 <- pars[2L]; pars_lin <- pars[3L]
  grd1 <- data.frame(K1=pars_K1/pars_lin,
                     K2=pars_K2)
  fit <- nls2(mod,
              start = grd1,
              algorithm = "brute-force",
              control = list(maxiter=100),
              lower=c(0,-10),upper=c(100,100))
  print(summary(fit))
  error <- function(y){
    pv  <- (predict(fit))
    rmse <- (rmse(y,predict(fit)))
    mae  <- (mae(y,predict(fit)))
    mse <- (mse(y,predict(fit)))
    rae <- (rae(y,predict(fit)))
    PAIC <- AIC(fit)
    PBIC <- BIC(fit)
    SE <- (sqrt(sum(predict(fit)-x)^2)/(n-2))
    colnames(y) <- rownames(y) <-colnames(y)
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
       main = "Hill-Deboer Isotherm Non-Linear Plot", xlab="theta", ylab= "Ce")
  lines(smooth.spline(x,predict(fit)),col="black")
  rsqq <- lm(Ce~predict(fit))
  print(summary(rsqq))
}

