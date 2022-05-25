#' @title Elovich Isotherm  Non-Linear Analysis
#' @name elovichanalysis
#' @description Elovich isotherm model is based on kinetic principle which
#' assumes that the adsorption sites would exponentially increase with chemical
#' reactions responsible for adsorption. It is suited for describing the behavior
#'  of adsorption concurring with the nature of chemisorption.
#' @param Ce the numerical value for equilibrium concentration
#' @param Qe the numerical value for adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Elovich isotherm, and
#' model error analysis
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples elovichanalysis(Ce,Qe)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Zeldowitsch, J. (1934). "Uber Den Mechanismus der Katalytischen
#' Oxidation Von CO a MnO2," URSS, Acta Physiochim, Vol. 1, No. 2, 1934, pp. 364-449.
#' @references Foo, K. Y., and Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export

# Building the Elovich isotherm nonlinear form

elovichanalysis <- function(Ce,Qe){

  x <- Qe
  y <- Ce
  data <- data.frame(x, y)

# Elovich isotherm nonlinear equation
  fit1 <- y ~ x / (QmE * KE * exp(-x/QmE))

# Setting of starting values
  start1 <- data.frame(KE = c(0, 100), QmE = c(0,100))

# Fitting of Elovich isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1,  data=data,
               control = nls.control(maxiter = 50, warnOnly = TRUE),
               algorithm = "port")

  print("Elovich Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Elovich iostherm model

  errors <- function(y) {
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
         "Standard Error for the Regression S" = SE )
  }
  a <- errors(y)
  print(a)

  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))

  # Graphical representation of the Elovich isotherm model

  ### Predicted parameter values
  parsElovich <- as.vector(coefficients(fit2))
  pars_KE <- parsElovich[1L];
  pars_QmE <- parsElovich[2L]

  rhs <- function (x){(x /(pars_QmE * pars_KE * exp(-x/pars_QmE)))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Qe",
         y = "Ce",
         title = "Elovich Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) + ggplot2::coord_flip()
}


