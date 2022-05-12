#' @title Dubinin-Radushkevich Isotherm  Non-Linear Analysis
#' @name dubininradushkevichanalysis
#' @description Dubinin-Radushkevich isotherm model is being utilized to define
#' adsorption energy mechanisms with Gaussian distribution onto heterogeneous surfaces.
#' Specifically, this model works well with an intermediate range of adsorbate
#' concentrations because it shows abnormal asymptotic behavior and is unable to
#' forecast Henry's Law at low pressure.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @param Temp temperature
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Dubinin-Radushkevich isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples Temp <- 298
#' @examples dubininradushkevichanalysis(Ce, Qe, Temp)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Dubinin, M.M. and Radushkevich, L.V. (1947) The Equation of the
#' Characteristic Curve of Activated Charcoal.
#' Proceedings of the Academy of Sciences, Physical Chemistry Section, 55, 331.
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export
#'



# Building the Dubinin-Radushkevich isotherm nonlinear form
dubininradushkevichanalysis <- function(Ce, Qe, Temp){

  x <- Ce
  y <- Qe
  R <- 8.314
  t <- Temp
  epsilon <- R*t*log(1+(1/x))
  data <- data.frame(x, y)

  ### Dubinin-Radushkevich isotherm nonlinear equation
  fit1 <- y ~ qs*exp(-K*(epsilon^2))

  ### Setting of starting values
  start1 <- list(qs = 1, K = 1e-7)

  ### Fitting of the Dubinin-Radushkevich isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1,  data=data,
               control = nls.control(maxiter = 100, warnOnly = TRUE),
               algorithm = "port")

  print("Dubinin-Radushkevich Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

  # Error analysis of the Dubinin-Radushkevich isotherm model

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


  # Graphical representation of the Dubinin-Radushkevich isotherm model

  ### Predicted parameter values
  parsdubR <- as.vector(coefficients(fit2))
  pars_qs <- parsdubR[1L];
  pars_K <- parsdubR[2L]

  rhs <- function(x){(pars_qs*exp(-pars_K*(R*t*log(1+(1/x)))^2))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Dubinin-Radushkevich Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}





