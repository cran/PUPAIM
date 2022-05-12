#' @title Radke-Prausnitz Isotherm Nonlinear Analysis
#' @name radkeprausnitzanalysis
#' @description The Radke-Prausnitz isotherm model has several important properties which provides a good fit over a wide range of adsorbate concentrations but more preferred in most adsorption systems at low adsorbate concentration.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Radke-Prausnitz isotherm, and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples radkeprausnitzanalysis(Ce,Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Radke, C. J. and Prausnitz, J. M. (1972) <doi:10.1021/i160044a003> Adsorption of organic solutions from dilute aqueous solution on activated carbon, Ind. Eng. Chem. Fund. 11 (1972) 445-451.
#' @export
#'

# Building the Radke-Parusnitz isotherm nonlinear form
radkeprausnitzanalysis <- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# Radke-Parusnitz isotherm nonlinear equation
  fit1 <- y ~ (Qmax*Krp*x)/(1+(Krp*x))^Mrp

# Setting of starting values
  N <- nrow(na.omit(data))
  start1 <- data.frame(Krp = seq(0, 100, length.out = N),
                       Qmax = seq(-1, 10, length.out = N),
                       Mrp = seq(-1, 1, length.out = N))

# Fitting of the Radke-Prausnitz isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 100 , warnOnly = TRUE),
                 algorithm = "port")

  print("Radke-Prausnitz Non-linear Analysis")
  print(summary(fit2))

  AIC <- AIC(fit2)
  print("Aikake Information Criterion")
  print(AIC)

  BIC <- BIC(fit2)
  print("Bayesian Information Criteron")
  print(BIC)

# Error analysis of the Radke-Prausnitz isotherm model

  errors <- function(y) {
    rmse <-Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
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

# Graphical representation of the Radke-Praustnitz isotherm model

  ### Predicted parameter values
  parsradkeP <- as.vector(coefficients(fit2))
  pars_Krp <- parsradkeP[1L];
  pars_Qmax <- parsradkeP[2L];
  pars_Mrp <- parsradkeP[3L]

  rhs <- function(x){((pars_Qmax*pars_Krp*x)/(1+(pars_Krp*x))^pars_Mrp)}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Radke-Praustnitz Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
