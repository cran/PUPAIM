#' @title Fritz-Schlunder Four Parameter Isotherm Non-Linear Analysis
#' @name FS4analysis
#' @description An empirical equation of Langmuir-Freundlich isotherm which
#' can fit a wide range of experimental results because of the large number of
#' coefficients in the isotherm.
#' @param Ce the numerical value for equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Fritz-Schlunder Four
#' Parameter isotherm, and model error analysis
#' @examples
#' \dontrun{
#'  Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600,0.63607, 0.80435, 1.10327, 1.58223)
#'  Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299,0.15379, 0.15735, 0.15735, 0.16607)
#'  FS4analysis(Ce,Qe)}
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Fritz, W., and Schluender, E. U. (1974) <doi:10.1016/0009-2509(74)80128-4> Simultaneous adsorption
#' equilibria of organic solutes in dilute aqueous solutions on activated carbon.
#' Chemical Engineering Science, 29(5), 1279-1282.
#' @export
#'

# Building the Fritz-Schlunder Four Parameter isotherm Nonlinear Form
FS4analysis <- function(Ce,Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# Fritz-Schlunder Four Parameter Nonlinear Equation

  fit1 <- y ~ (C*x^a)/(1+D*x^b)

# Setting of Starting Values
  start <- data.frame(C = c(0, 1000),
                      D = c(0, 1000),
                      a = c(0, 1),
                      b = c(0, 1))

# Fitting of the Fritz-Schlunder Four Parameter isotherm via nls2
  fit2 <- nls2::nls2(fit1, start = start, data = data,
               control = nls.control(maxiter = 10000 , warnOnly = TRUE),
               algorithm = "plinear-random")

  pars <- as.vector(coefficients(fit2))
  pars_C <- pars[1L];
  pars_D <- pars[2L];
  pars_a <- pars[3L];
  pars_b <- pars[4L];
  pars_lin <- pars[5L]

  start1 <- data.frame(C = pars_C*pars_lin,
                       D = pars_D,
                       a = pars_a,
                       b = pars_b)

  fit3 <- nls2::nls2(fit1,
               start = start1, data = data, control=list(maxiter= 100, warnOnly = T),
               algorithm = "default")


  print("Fritz-Schlunder Four Parameter Isotherm Nonlinear Analysis")
  print(summary(fit3))

  print("Akaike Information Criterion")
  print(AIC(fit3))

  print("Bayesian Information Criterion")
  print(BIC(fit3))

# Error Analysis of the Fritz-Schlunder Four Parameter Isotherm Model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit3))
    mae <- Metrics::mae(y, predict(fit3))
    mse <- Metrics::mse(y, predict(fit3))
    rae <- Metrics::rae(y, predict(fit3))
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

  rsqq <- lm(Qe~predict(fit3))
  print(summary(rsqq))

# Graphical representation of the Fritz-Schlunder Four Parameter isotherm model

  ### Predicted parameter values
  parsfritzV <- as.vector(coefficients(fit3))
  pars_C <- parsfritzV[1L];
  pars_D <- parsfritzV[2L];
  pars_a <- parsfritzV[3L];
  pars_b <- parsfritzV[4L];

  rhs <- function(x){((pars_C*x^pars_a)/(1+pars_D*x^pars_b))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Fritz-Schlunder (IV) Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
