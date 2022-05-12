#' @title Brunauer-Emett-Teller (BET) Isotherm Non-Linear Analysis
#' @name BETanalysis
#' @description BET was particularly formulated to describe the multilayer adsorption
#' process in gas systems, but can also be employed to an aqueous solution that relates
#' the binding between layers because of the molecular charge among them.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for BET isotherm, and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples BETanalysis(Ce,Qe)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Brunauer, S., Emmett, P.H. and Teller, E. (1938) <doi:10.1021/ja01269a023> Adsorption of Gases in Multimolecular Layers.
#' Journal of the American Chemical Society, 60, 309-319.
#' @export
#'

# Building the BET nonlinear form
BETanalysis <- function(Ce,Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# BET isotherm nonlinear equation
  fit1 <- y ~ (CBET *x)/((Cs-x)*(1+((CBET-1)*(x/Cs)))) ### Qmax is conditionally linear

# Setting of starting values
  N <- nrow(na.omit(data))
  start1 <- data.frame(CBET = seq(-50, 400, length.out = 50),
                       Cs = seq(-400, 50, length.out = 50))

# Fitting of the BET isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1,  data=data,
                control = nls.control(maxiter = 45 , warnOnly = TRUE),
                algorithm = "plinear-random")

  print("Brunauer-Emett-Teller (BET) Isotherm Non-Linear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

#Error analysis of the BET Isotherm model

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

  # Graphical representation of the BET isotherm model

  ### Predicted parameter values
  parsBET <- as.vector(coefficients(fit2))
  pars_CBET <- parsBET[1L];
  pars_Cs <- parsBET[2L];
  pars_Qmax <- parsBET[3L]

  rhs <- function (x){(pars_Qmax*pars_CBET*x)/((pars_Cs-x)*(1+((pars_CBET-1)*(x/pars_Cs))))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "BET Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

