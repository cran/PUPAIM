#' @title Aranovich Isotherm Non-Linear Analysis
#' @name aranovichanalysis
#' @description The Aranovich isotherm (Aranovich, 1992) is a three-parameter
#' isotherm model that is a modified version of the BET isotherm. This isotherm
#' model is theoretically corrected by polymolecular adsorption isotherm and is
#' applicable to modeling adsorption with a wide range concentration of the
#' adsorbate molecules.
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
#' @examples aranovichanalysis(Ce,Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Aranovich, G. L. (1992) <doi:10.1021/la00038a071> The Theory of Polymolecular Adsorption.
#' Langmuir, 8(2), 736-739.
#' @export
#'

# Building the Aranovich isotherm nonlinear form
aranovichanalysis <- function(Ce,Qe) {

  x <- Ce
  y <- Qe
  data <- data.frame(Ce,Qe)

# Aranovich isotherm nonlinear equation
  fit1 <- y ~ (Qmax*CA*(x/CsA))/(sqrt(1-(x/CsA))*(1+(CA*(x/CsA))))

# Setting of starting values
  start1 <- data.frame(Qmax = c(1,10) , CA = c(1, 100), CsA = c(10, 100))

# Fitting of the Aranovich isotherm via nls2


  suppressWarnings(fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 100 , warnOnly = TRUE),
                 algorithm = "default"))

  print("Aranovich Isotherm Non-Linear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

#Error analysis of the Aranovich Isotherm model

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
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))
  
# Graphical representation of the Aranovich isotherm model

  ### Predicted parameter values
  parsara <- as.vector(coefficients(fit2))
  pars_Qmax <- parsara[1L];
  pars_ca <- parsara[2L];
  pars_csa <- parsara[3L]

  rhs <- function(x){(pars_Qmax*pars_ca*(x/pars_csa))/
      (sqrt(1-(x/pars_csa))*(1+(pars_ca*(x/pars_csa))))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Aranovich Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}
