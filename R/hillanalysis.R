#' @title Hill Isotherm Non-Linear Analysis
#' @name hillanalysis
#' @description Hill isotherm model shows the connection of different species
#' being adsorbed on to the homogeneous surfaces. This isotherm model supposes
#' that adsorption is a cooperative phenomenon which means the adsorbates having
#' the capability to bind at one specific site on the adsorbent affecting other
#' binding sites on the same adsorbent
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Hill isotherm, and model
#' error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples hillanalysis(Ce,Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Hill, T. L. (1946) <doi:10.1063/1.1724129> "Statistical mechanics of multimolecular
#' adsorption II. Localized and mobile adsorption and absorption," The Journal
#' of Chemical Physics, vol. 14, no. 7, pp. 441-453.
#' @export

# Building the Hill isotherm nonlinear form
hillanalysis <- function(Ce,Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

  ### Hill isotherm nonlinear equation
  fit1 <- (y ~ ((qh*x^nh)/(Kd+x^nh)))

  ### Setting of starting values
  start1 <- list(qh = 1, nh = 1, Kd = 1)

  ### Fitting of the Hill isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 100, warnOnly = TRUE),
                 algorithm = "port")

  print("Hill Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

  # Error analysis of the Hill isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(x, predict(fit2))
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
  
  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))
  
# Graphical representation of the Hill isotherm model

  ### Predicted parameter values
  parshill <- as.vector(coefficients(fit2))
  pars_qh <- parshill[1L];
  pars_nh <- parshill[2L];
  pars_Kd <- parshill[3L]

  rhs <- function(x){((pars_qh*x^pars_nh)/(pars_Kd+x^pars_nh))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
       y = "Qe",
       title = "Hill Isotherm Nonlinear Model",
       caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
