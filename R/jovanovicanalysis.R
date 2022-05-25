#' @title Jovanovic Isotherm Non-Linear Analysis
#' @name jovanovicanalysis
#' @description The Jovanovic isotherm model was built upon the assumptions
#' based on the Langmuir isotherm model with few possible inclusions of
#' mechanical contact among the desorbing and adsorbing molecules. The adjustment
#' of the adsorption surface from this model made the equation less effective
#' in the physical adsorption but can be applied to adsorption with both mobile
#' and localized monolayer without lateral interaction. Moreover, the equation
#' of the Jovanovic isotherm model is able to reach the limit of saturation when
#' there is high concentration, while it reduces to Henry's Law at low concentration.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorpted capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Jovanovic isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples jovanovicanalysis(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references:  Saadi, R., Saadi, Z., Fazaeli, R., Fard, N. E. (2015) <DOI: 10.1007/s11814-015-0053-7>
#' Monolayer and multilayer adsorption isotherm models for sorption from aqueous media.
#' Korean J. Chem. Eng., 32(5), 787-799 (2015)
#' @references:  Vargas, A., Cazetta, A., Kunita, M., Silva, T., Almeida V. (2011) <DOI:10.1016/j.cej.2011.01.067>
#' Adsorption of methylene blue on activated carbon produced from Flamboyant pods
#' (Delonix regia): Study of adsorption isotherms and kinetic models. Chemical
#' Engineering Journal 168 (2011) 722-730
#' @export
#'

# Building the Jovanovic isotherm nonlinear form
jovanovicanalysis<- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data<- data.frame(x, y)

  # Jovanovic isotherm nonlinear equation
  fit1 <- y ~ (Qmax*(1 - exp(-Kf*x)))

  # Setting of starting values
  start1 <- list(Qmax= 1, Kf= 1)

  # Fitting of the Jovanovic isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1,  data=data,
                 control = nls.control(maxiter = 50, warnOnly = TRUE),
                 algorithm = "port")

  print("Jovanovic Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Aikake Information Criterion")
  print(AIC(fit2))

  BIC <- BIC(fit2)
  print("Bayesian Information Criterion")
  print(BIC(fit2))

  # Error analysis of the Jovanovic isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit2))
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

  # Graphical representation of the Jovanovic isotherm model

  ### Predicted parameter values
  parsjova <- as.vector(coefficients(fit2))
  pars_Qmax <- parsjova[1L];
  pars_Kf <- parsjova[2L];

  rhs <- function(x) {((pars_Qmax*(1 - exp(-pars_Kf*x))))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Jovanovic Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}

