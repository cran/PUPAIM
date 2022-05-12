#' @title Baudu Isotherm Non-Linear Analysis
#' @name bauduanalysis
#' @description Baudu is a reduced form of Langmuir isotherm since it was observed
#' that the estimation of Langmuir coefficients b and qm by tangent measurements
#' at different equilibrium constants are not constants in the broad concentration
#'  range. This can be used if the ranges are (1+x+y) <1 and (1+x) <1.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Baudu isotherm, and model error
#' analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples bauduanalysis(Ce,Qe)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Baudu, M. (1990). Etude des interactions solutes-fibres de charbon actif:
#' applications et regeneration (Doctoral dissertation, Rennes 1).
#' from https://www.theses.fr/1990REN10039
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export
#'

# Building the Baudu isotherm nonlinear form
bauduanalysis <- function(Ce,Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# Baudu isotherm nonlinear equation
  fit1 <- y ~ ((Qmax*bo*x)^(1+a+b))/(1+bo*x^(1+a))

# Setting of starting values
  start1 <- data.frame(Qmax= c(0,100), bo= c(0,100), a= c(0, 10), b= c(-10,1))

# Fitting of Baudu isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data = data,
                 control= nls.control(maxiter=45, warnOnly=TRUE),
                 algorithm= "port")

  print("Baudu Isotherm Non-linear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print((fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))


# Error analysis of the Baudu isotherm model

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
         "Root Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Baudu isotherm model

  ### Predicted parameter values
  parsbaudu <- as.vector(coefficients(fit2))
  pars_Qmax <- parsbaudu[1L];
  pars_bo <- parsbaudu[2L];
  pars_a <- parsbaudu[3L];
  pars_b <- parsbaudu[4L];

  rhs <- function(x){((pars_Qmax*pars_bo*x)^(1+pars_a+pars_b))/(1+pars_bo*x^(1+pars_a))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Baudu Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
