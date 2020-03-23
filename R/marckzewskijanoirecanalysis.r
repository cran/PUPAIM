#' @title Marckzewski-Jaroniec Isotherm Analysis Non-Linear Form
#' @description The Marczewski-Jaroniec Isotherm model is the resemblance of Langmuir Isotherm model. It is developed on the basis of the supposition of local Langmuir isotherm and adsorption energies distribution in the active sites on adsorbent (Parker, 1995; Sivarajasekar & Baskar, 2014). This equation comprises all isotherm equations being an extension of simple Langmuir Isotherm to single solute adsorption on heterogeneous solids (Marczewski & Jaroniec, 1983).
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the absorbed capacity
#' @importFrom graphics "abline" "plot"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return The non-linear regression and the parameters for Marckzewski Jaroniec Isotherm Analysis
#' @examples mjanalysis(moringa$Ce,moringa$Qe)
#' @author Brent Mark A. Daniel
#' @author Charlestone E. Antatico
#' @author C.C. Deocaris
#' @references Marczewski, A.W., Jaroniec, M. (1983). A new isotherm equation for single-solute adsorption from dilute solutions on energetically heterogeneous solids. Monatsh Chem 114, 711-715. doi: 10.1007/BF01134184
#' @references Parker Jr, G.R. (1995). Optimum isotherm equation and thermodynamic interpretation for aqueous 1,1,2- trichloroethene adsorption isotherms on three adsorbents. Adsorption, 1 (2): 113-132. doi:10.1007/BF00705000
#' @references Sivarajasekar, N., Baskar, R. (2014). Adsorption of basic red 9 onto activated carbon derived from immature cotton seeds: Isotherm studies and error analysis. Desalination and Water Treatment, 52: 1-23. doi:10.1080/19443994.2013.834518
#' @export
mjanalysis<- function (Ce, Qe)
{
  .data <- data.frame(Ce,Qe)
  x <- Ce
  y <- Qe
  mod252 <- (y ~ ((q*((K * x)^n)) /(1+((K * x)^n)))^(m/n))
  K_max <- ((mean(y)/mean(x))+1)
  N <- nrow(na.omit(.data))
  values <- data.frame(q = seq(0, 1000, length.out = N), K = seq(0, K_max, length.out = N), m = seq(0, 1, length.out = N), n = seq(0, 1, length.out = N))
  set.seed(511)
  fit253 <- nls2(mod252, .data, start = values, algorithm  = "plinear-random", control = nls.control(maxiter = 1000, tol = 1e-05, warnOnly = T),
                low = c(0, 0, 0, 0), upper = c(1000, K_max, 1, 1))
  pars <- as.vector(coefficients(fit253))
  pars_q <- pars[1L]; pars_K <- pars[2L]; pars_m <- pars[3L]; pars_n <- pars[4L]
  print("Non-Linear Analysis For Marckzewski Jaroniec Isotherm")
  print(summary(fit253))
  error <- function(y){
    pv  <- (predict(fit253))
    rmse<- (rmse(y,predict(fit253)))
    mae <- (mae(y,predict(fit253)))
    mse <- (mse(y,predict(fit253)))
    rae <- (rae(y,predict(fit253)))
    PAIC <- AIC(fit253)
    PBIC <- BIC(fit253)
    SE <-(sqrt(sum(predict(fit253)-x)^2)/(N-2))
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
    non_linear_plot <- function (Ce, Qe)
  {
    plot(x, y, main = "Marckzewski-Jaroniec Isotherm Non-Linear Plot", xlab = "Ce", ylab = "Qe")

    lines(smooth.spline(x,predict(fit253)), col = "black")
  }
  non_linear_plot (Ce,Qe)
  rsqq <- lm(Qe~predict(fit253))
  print(summary(rsqq))

}
