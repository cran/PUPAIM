#'@title Langmuir Isotherm Nonlinear Analysis
#'@name temkinanalysis
#'@description The Langmuir isotherm is described to be the most useful and
#' simplest isotherm for both chemical adsorption and physical adsorption. It
#' assumes that there is uniform adsorption energy onto the monolayer surface
#' and that there would be no interaction between the adsorbate and the surface.
#'@param Ce the numerical value for the equilibrium capacity
#'@param Qe the numerical value for the adsorbed capacity
#'@import nls2
#'@import Metrics
#'@import stats
#'@import ggplot2
#'@return the nonlinear regression, parameters for Langmuir isotherm, and model
#'error analysis
#'@examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#'@examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#'@examples langmuiranalysis(Ce,Qe)
#'@author Keith T. Ostan
#'@author Chester C. Deocaris
#'@references Langmuir, I. (1918) <doi:10.1021/ja01269a066> The adsorption of gases on plane surfaces of
#'glass, mics and platinum. Journal of the American Chemical Society, 1361-1403.
#'@export
#'

# Building the Langmuir isotherm nonlinear form
langmuiranalysis <- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# Langmuir isotherm nonlinear equation
  fit1 <- y ~ (Qmax*Kl*x)/(1+(Kl*x))

# Setting of starting values
  N <- nrow(na.omit(data))
  start1<- list(Qmax = 1, Kl= 1)

# Fitting of the Langmuir isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter=50, warnOnly = TRUE),
                 algorithm = "port")

  print("Langmuir Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Langmuir isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Root Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Langmuir isotherm model

  ### Predicted parameter values
  parslang <- as.vector(coefficients(fit2))
  pars_Kl <- parslang[2L];
  pars_Qmax <- parslang[1L];

  rhs <- function(x){((pars_Qmax*pars_Kl*(x))/(1+(pars_Kl*x)))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Langmuir Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}

