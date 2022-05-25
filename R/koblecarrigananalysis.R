#' @title Koble-Carrigan Isotherm Nonlinear Analysis
#' @name koblecarrigananalysis
#' @description It is three-parameter isotherm model equation that incorporates
#' both Freundlich and Langmuir isotherms for representing equilibrium adsorption
#' data. Koble-Corrigan isotherm model appeared to have advantages over both the
#' Langmuir and Freundlich equations in that it expresses adsorption data over
#' very wide ranges of pressures and temperatures.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Koble-Carrigan isotherm, and
#' model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples koblecarrigananalysis(Ce, Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Corrigan, T. E., and Koble, R. A.(1952) <doi:10.1021/ie50506a049> Adsorption isotherms for
#' pure hydrocarbons Ind. Eng. Chem. 44 383-387.
#' @export

# Building the Koble-Corrigan isotherm nonlinear form
koblecarrigananalysis <- function(Ce,Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# Koble-Corrigan isotherm nonlinear equation
  fit1 <- y ~ (Ak*(x^p))/(1 + Bk*(x^p))

# Setting of starting values
  start1 <- list(Ak = 1, Bk = 1, p = 1)

# Fitting of the Koble-Corrigan isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter =50 , warnOnly = TRUE),
                 algorithm = "port")

  print("Koble-Carrigan Isotherm Nonlinear Analysis")
  print(summary(fit2))

  AIC <- AIC(fit2)
  print("Aikake Information Criterion")
  print(AIC)

  BIC <- BIC(fit2)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Koble-Corrigan isotherm model

  errors <- function(y) {
  rmse <- Metrics::rmse(y, predict(fit2))
  mae <- Metrics::mae(y, predict(fit2))
  mse <- Metrics::mse(y, predict(fit2))
  rae <- Metrics::rae(y, predict(fit2))
  N <- nrow(na.omit(data))
  SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
  colnames(y) <- rownames(y) <- colnames(y)
  list("Relative Mean squared Error" = rmse,
       "Mean Absolute Error" = mae,
       "Mean Squared Error" = mse,
       "Relative Absolute Error" = rae,
       "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)
  
  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))

# Graphical representation of the Koble-Corrigan isotherm model

  ### Predicted parameter values
  parskobleC <- as.vector(coefficients(fit2))
  pars_Ak <- parskobleC[1L];
  pars_Bk <- parskobleC[2L];
  pars_p <- parskobleC[3L]

  rhs <- function(x){((pars_Ak*(x^pars_p))/(1 + pars_Bk*(x^pars_p)))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Koble-Corrigan Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
