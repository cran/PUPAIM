#'@title Toth Isotherm  Nonlinear Analysis
#'@name temkinanalysis
#'@description Another empirical modification of the Langmuir equation with the
#'aim of reducing the error between experimental data and predicted value of
#'equilibrium data.
#'@param Ce the numerical value for the equilibrium capacity
#'@param Qe the numerical value for the fractional coverage
#'@import nls2
#'@import Metrics
#'@import stats
#'@import ggplot2
#'@return the nonlinear regression, parameters for Toth isotherm, and
#'model error analysis
#'@examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#'@examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#'@examples tothanalysis(Ce,Qe)
#'@author Keith T. Ostan
#'@author Chester C. Deocaris
#'@references Toth, J. (1971). State equations of the solid gas interface layer.
#'Acta Chem. Acad. Hung. 69:311-317
#'@export
#'

# Building the Toth isotherm nonlinear form
tothanalysis <- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x,y)

# Toth isotherm nonlinear equation
  fit1 <- y ~ (x)/(At+x)^(1/Nt) ##Kt is conditionally linear

# Setting of starting values
  start1 <- data.frame(At = c(1, 100), Nt = c(0, 1))

# Fitting of the Toth isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 100, warnOnly = TRUE),
                 algorithm = "plinear-random")


  print("Toth Isotherm Nonlinear Analysis")
  print(summary(fit2))

  AIC <- AIC(fit2)
  print("Aikake Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))


# Error analysis of the Sips isotherm model

  errors <- function (y) {
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    N <- nrow(na.omit(data))
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

# Graphical representation of the Toth isotherm model

  ### Predicted parameter values
  parstoth <- as.vector(coefficients(fit2))
  pars_At <- parstoth[1L];
  pars_Nt <- parstoth[2L];
  pars_Kt <- parstoth[3L]

  rhs <- function(x){((pars_Kt*(x))/(pars_At+x)^(1/pars_Nt))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Toth Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
