#' @title Volmer Isotherm Non-Linear Analysis
#' @name volmeranalysis
#' @description The Volmer isotherm describes a distribution of monolayer
#' adsorption processes. This theoretical model has the assumption in which
#' the adsorbate molecules can move toward the surfaces of adsorbents, and
#' the interactions that can be formed between the adsorbates are negligible.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Aranovich isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples volmeranalysis(Qe,Ce)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Volmer, M. (1925) <doi:10.1515/zpch-1925-11519> Thermodynamische folgerungen aus der
#' zustandsgleichung fur adsorbierte stoffe. Z. Phys. Chem. 115, 253-261.
#' @export
#'

# Building the Volmer isotherm nonlinear form
volmeranalysis <- function(Qe,Ce) {

  x <- Qe
  y <- Ce
  data <- data.frame(Qe,Ce)

  # Volmer isotherm nonlinear equation
  fit1 <- y ~ ((1/bV)*(x/(Qmax - x))*exp((x/(Qmax - x))))

  # Setting of starting values
  start1 <- data.frame(Qmax = c(1, 1000), bV = c(1, 100))

  # Fitting of the Volmer isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 100, warnOnly = TRUE),
                 algorithm = "port")

  print("Volmer Isotherm Non-Linear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))


  #Error analysis of the Volmer Isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    N <- nrow(na.omit(data))
    SE <- SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Root Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  # Graphical representation of the Volmer isotherm model

  ### Predicted parameter values
  parsvol <- as.vector(coefficients(fit2))
  pars_Qmax <- parsvol[1L];
  pars_bV <- parsvol[2L];

  rhs <- function(x){((1/pars_bV)*(x/(pars_Qmax - x))*exp((x/(pars_Qmax - x))))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Qe",
         y = "Ce",
         title = "Volmer Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}
