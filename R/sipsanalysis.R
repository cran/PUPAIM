#'@title Sips Isotherm Nonlinear Analysis
#'@name sipsanalysis
#'@description It is the most applicable to use in the monolayer adsorption
#'isotherm model amongst the three-parameter isotherm models and is also valid
#'for the prediction of heterogeneous adsorption systems as well as localized
#'adsorption with no interactions occurring between adsorbates.
#'@param Ce the numerical value for the equilibrium capacity
#'@param Qe the numerical value for the adsorbed capacity
#'@import nls2
#'@import Metrics
#'@import stats
#'@import ggplot2
#'@return the nonlinear regression, parameters for Sips isotherm, and model
#'error analysis
#'@examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#'@examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#'@examples sipsanalysis(Ce,Qe)
#'@author Keith T. Ostan
#'@author Chester C. Deocaris
#'@references Sips, R. (1948) <doi:10.1063/1.1746922> On the structure of a catalyst surface.
#'The Journal of Chemical Physics, 16(5), 490-495.
#'@export
#'

# Building the Sips isotherm nonlinear form
sipsanalysis <- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data<- data.frame(x, y)

# Sips isotherm nonlinear equation
  fit1 <- y ~ (Ks*(x^n))/(1 + (As*(x^n)))

# Setting of starting values
  start1 <- list(As = 1, Ks = 1, n = 1)

# Fitting of the Sips isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 50, warnOnly = TRUE),
                 algorithm = "port")

  print("Sips Isotherm Nonlinear Analysis")
  print(summary(fit2))

  AIC <- AIC(fit2)
  print("Aikake Information Criterion")
  print(AIC)

  BIC <- BIC(fit2)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Sips isotherm model

errors <- function(y){
    rmse <- Metrics::rmse(Qe, predict(fit2))
    mae <- Metrics::mae(Qe, predict(fit2))
    mse <- Metrics::mse(Qe, predict(fit2))
    rae <- Metrics::rae(Qe, predict(fit2))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae)
  }
  s <- errors(y)
  print(s)

# Graphical representation of the Sips isotherm model

  ### Predicted parameter values
  parssips <- as.vector(coefficients(fit2))
  pars_As <- parssips[1L];
  pars_Ks <- parssips[2L];
  pars_n <- parssips[3L]

  rhs <- function(x){((pars_Ks*(x^pars_n))/(1 + (pars_As*(x^pars_n))))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Sips Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
