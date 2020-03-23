#' @title Halsey Isotherm Analysis Non-Linear Form
#' @description used to evaluate multilayer adsorption at a relatively large distance from the surface
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the nonlinear regression and the parameters for the Halsey isotherm analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples halseyanalysis(Ce, Qe)
#' @author Aries N. Bunag
#' @author C.C. Deocaris
#' @references Sousa Neto, V.O., Oliveira, A. G., Teixeira, R.N.P., Silva, M.A.A.,Freire, P.T.C., Keukeleire, D.D., & Nascimento, R.S.(2011). USE OF COCONUT BAGASSE AS
#' ALTERNATIVE ADSORBENT FOR SEPARATION OF COPPER(III) IONS FROM AQUEOUS SOLUTIONS: ISOTHERMS, KINETIC AND THERMODYNAMIC STUDIES. Retrived February 17, 2020,
#' from https://bioresources.cnr.ncsu.edu/BioRes_06/BioRes_06_3_3376_Neto_OTSFKN_Coconut_Bagasse_Ads_Cu2_Water_Kinet_Thermo_1822.pdf
#' @references Imran, M., Naseem, Khalida, Mirza, Latif, M., & Madeeha. (2018, December 1). Evaluation of Saccharum bengalense as a Non-Conventional Biomaterial for
#' Biosorption of Mn(II) Ions from Aqueous Solutions. Retrieved February 17.2020, from http://www.ijcce.ac.ir/article_29361.html
#' @export
halseyanalysis <- function(Ce,Qe)
{
  x <- Ce
  y <- Qe
  data <- data.frame(x, y)
  mod224 <- Qe ~ (Kh/Ce)^(1/nh)
  N <- nrow(na.omit(data))
  values <- data.frame(nh = seq(-10, 10, length.out = N), Kh = seq(0, 1000, length.out = N))
  set.seed(511)
  suppressWarnings(fit225 <- nls2(mod224, data = data, start= values, control = nls.control(maxiter = 1000 , warnOnly = TRUE), algorithm = "port"))
  print("NLS2 Analysis for Halsey Isotherm Analysis")
  print(summary(fit225))
  error <- function(y){
    pv  <- (predict(fit225))
    rmse<- (rmse(y,predict(fit225)))
    mae <- (mae(y,predict(fit225)))
    mse <- (mse(y,predict(fit225)))
    rae <- (rae(y,predict(fit225)))
    PAIC <- AIC(fit225)
    PBIC <- BIC(fit225)
    SE <-(sqrt(sum(predict(fit225)-x)^2)/(N-2))
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
    plot(x, y, main = "Halsey Analysis", xlab = "Ce",
         ylab = "Qe")
    lines(x, predict(fit225), col = "black")
    rsqq <- lm(Qe~predict(fit225))
    print(summary(rsqq))
  }
