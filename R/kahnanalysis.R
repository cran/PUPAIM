#' @title Kahn Isotherm Non-Linear Analysis
#' @name kahnanalysis
#' @description A generalized model recommended for pure solutions, in which
#' both extremes, Langmuir and Freundlich, can be represented. This isotherm was
#' developed to cater to both the single- and multi-component adsorption systems.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value the absorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Kahn isotherm, and
#' model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples kahnanalysis(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Khan, A. R., Al-Waheab, I. R., & Al-Haddad, A. (1996) <doi:10.1080/09593331708616356> A
#' generalized equation for adsorption isotherms for multi-component organic
#' pollutants in dilute aqueous solution. Environmental Technology (United Kingdom), 17(1), 13-23.
#' @export
#'

# Building the Kahn isotherm nonlinear model
kahnanalysis <- function(Ce,Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(Ce, Qe)

# Kahn isotherm nonlinear equation
 fit1 <- y ~ (Qmax*bk*x)/((1+bk*x)^ak)

# Setting of starting values
 N <- nrow(na.omit(data))
 start1 <- list(Qmax = seq(1, 1000, length.out = N),
                bk = seq(1, 10, length.out = N),
                ak = seq(1, 10, length.out = N))

# Fitting of the Kahn isotherm via nls2

 fit2 <- nls2::nls2(fit1, start = start1, data=data,
           control = nls.control(maxiter = 1000 , warnOnly = TRUE),
           algorithm = "port")

 print("Kahn Isotherm Nonlinear Analysis")
 print(summary(fit2))

 print("Akaike Information Criterion")
 print(AIC(fit2))

 print("Bayesian Information Criterion")
 print(BIC(fit2))

  # Error analysis of the Kahn isotherm model

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

# Graphical representation of the Kahn isotherm model

  ### Predicted parameter values
  parskahn <- as.vector(coefficients(fit2))
  pars_Qmax <- parskahn[1L];
  pars_bk <- parskahn[2L];
  pars_ak <- parskahn[3L]

  rhs <- function(x){((pars_Qmax*pars_bk*x)/((1+pars_bk*x)^pars_ak))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Kahn Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
