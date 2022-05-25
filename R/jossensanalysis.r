#' @title Jossens Isotherm Non-Linear Analysis
#' @name jossensanalysis
#' @description The Jossens isotherm model predicts a simple equation based on
#' the energy distribution of adsorbate-adsorbent interactions at adsorption sites.
#' This model assumes that the adsorbent has heterogeneous surface with respect
#' to the interactions it has with the adsorbate.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Jossens isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples jossensanalysis(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Jossens, L., Prausnitz, J. M., Fritz, W., Schlunder, E. U.,
#' and Myers, A. L. (1978) <doi:10.1016/0009-2509(78)85015-5> Thermodynamics of
#' multi-solute adsorption from dilute aqueous solutions.
#' Chemical Engineering Science, 33(8), 1097-1106.
#' @export
#'

# Building the Jossens isotherm nonlinear model
jossensanalysis <- function(Ce, Qe){

 x <- Qe
 y <- Ce
 data <- data.frame(x, y)

# Jossens isotherm nonlinear equation
 fit1 <- y ~ (x/H)*(exp(J*(x^p)))

# Setting of starting values
 start1 <- data.frame(H = c(1, 500), J = c(1, 100), p = c(0, 10))

# Fitting of the Jossens isotherm via nls2

 fit2 <- nls2::nls2(fit1, start = start1,  data=data,
                control = nls.control(maxiter = 100 , warnOnly = TRUE),
                algorithm = c("port"))

 print("Jossens Isotherm Nonlinear Analysis")
 print(summary(fit2))

 print("Akaike Information Criterion")
 print(AIC(fit2))

 print("Bayesian Information Criterion")
 print(BIC(fit2))

# Error analysis of the Jossens isotherm model

errors <- function(y) {
  rmse <- Metrics::rmse(y, predict(fit2))
  mae <- Metrics::mae(y, predict(fit2))
  mse <- Metrics::mse(y, predict(fit2))
  rae <- Metrics::rae(y, predict(fit2))
  N <- nrow(na.omit(data))
  SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
  colnames(y) <- rownames(y) <- colnames(y)
  list( "Root Mean Squared Error" = rmse,
        "Mean Absolute Error"=mae ,
        "Mean Squared Error"=mse ,
        "Relative Absolute Error"=rae ,
        "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))

# Graphical representation of the Jossens isotherm model

  ### Predicted parameter values
  parsjoss <- as.vector(coefficients(fit2))
  pars_H <- parsjoss[1L];
  pars_J <- parsjoss[2L];
  pars_p <- parsjoss[3L]

  rhs <- function(x){((x/pars_H)*(exp(pars_J*(x^pars_p))))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Qe",
         y = "Ce",
         title = "Jossens Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5)) + ggplot2::coord_flip()
}
