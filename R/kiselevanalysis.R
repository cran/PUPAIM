#' @title Kiselev Isotherm Non linear Analysis
#' @name kiselevanalysis
#' @description It is also known as localized monomolecular layer model and is
#' only valid for surface coverage theta > 0.68.
#' @param theta is the fractional surface coverage
#' @param Ce the numerical value for equilibrium capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Kiselev isotherm, and
#' model error analysis
#' @examples theta <- c(0.19729, 0.34870, 0.61475, 0.74324, 0.88544, 0.89007, 0.91067, 0.91067, 0.96114)
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples kiselevanalysis(Ce, theta)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Kiselev, A. V. (1958). "Vapor adsorption in the formation of
#' adsorbate molecule complexes on the surface," Kolloid Zhur, vol. 20, pp. 338-348.
#' @export

# Building the Kiselev isotherm nonlinear model
kiselevanalysis <- function(Ce, theta){

  x <- theta
  y <- Ce
  data <- data.frame(x, y)

  # Kiselev isotherm nonlinear equation
  fit1 <- y ~ (x/(Ki*(1-x) * (1 + Kn*x)))

  # Setting of starting values
  start1 <- data.frame(Ki = c(-100,1000), Kn = c(-100,1000))

  # Fitting of the Kiselev isotherm via nls2

  suppressWarnings(fit2 <-  nls2::nls2(fit1, start = start1,  data=data,
                                  control = nls.control(maxiter= 100, warnOnly = TRUE),
                                  algorithm = "port"))

  print("Kiselev Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Kiselev isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae (y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y)<- rownames(y) <- colnames(y)
    list("Root Mean Square Error"= rmse,
       "Mean Absolute Error"= mae,
       "Mean Squared Error"= mse,
       "Relative Absolute Error"= rae,
       "Standard Error for the Regression S" = SE)
    }
    a <- errors(y)
    print(a)
    
    rsqq <- lm(theta~predict(fit2))
    print(summary(rsqq))

# Graphical representation of the Kiselev isotherm model

  ### Predicted parameter values
  parskise <- as.vector(coefficients(fit2))
  pars_Ki <- parskise[1L];
  pars_Kn <- parskise[2L];

  rhs <- function(x){((x/(pars_Ki*(1-x) * (1 + pars_Kn*x))))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = expression(paste(theta)),
         y = "Ce",
         title = "Kiselev Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5)) + ggplot2::coord_flip()
}
