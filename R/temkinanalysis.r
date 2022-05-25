#' @title Temkin Isotherm Nonlinear Analysis
#' @name temkinanalysis
#' @description Temkin isotherm  is a monolayer adsorption isotherm model which
#' takes into account the effects that the indirect interaction amongst adsorbate
#' molecules could have on the adsorption process.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @param Temp temperature
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Temkin isotherm, and model
#' error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples Temp <- 298
#' @examples temkinanalysis(Ce, Qe, Temp)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Temkin, M.J., and Pyzhev, V. (1940). Kinetics of ammonia synthesis
#' on promoted iron catalyst. Acta Phys. Chim. USSR 12, 327-356.
#' @export
#'

# Building the Temkin isotherm nonlinear form
temkinanalysis <- function(Ce, Qe, Temp){

  x <- Ce
  y <- Qe
  t <- Temp
  R <- 8.314
  data<- data.frame(x, y)

# Temkin isotherm nonlinear equation
  fit1 <- y ~ ((R*t)/bT) * log(AT*x)

# Setting of starting values
  start1 <- list(AT= 1, bT= 1)

# Fitting of the Temkin isotherm via nls2


  fit2 <- nls2::nls2(fit1, start = start1, data=data,
              control = nls.control(maxiter = 50, warnOnly = TRUE),
              algorithm = "default")

  print("Temkin Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Aikake Information Criterion")
  print(AIC(fit2))

  BIC <- BIC(fit2)
  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Temkin isotherm model

errors <- function(y) {
  rmse <-Metrics::rmse(y, predict(fit2))
  mae <- Metrics::mae(y, predict(fit2))
  mse <- Metrics::mse(y, predict(fit2))
  rae <- Metrics::rae(y, predict(fit2))
  N <- nrow(na.omit(data))
  SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
  colnames(y) <- rownames(y) <- colnames(y)
  list("Relative Mean Squared Error" = rmse,
       "Mean Absolute Error" = mae,
       "Mean Squared Error" = mse,
       "Relative Absolute Error" = rae,
       "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)
  
  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))

# Graphical representation of the Temkin isotherm model

  ### Predicted parameter values
  parstem <- as.vector(coefficients(fit2))
  pars_AT <- parstem[1L];
  pars_bT <- parstem[2L];

  rhs <- function(x){(((R*t)/pars_bT) * log(pars_AT*x))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Temkin Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
