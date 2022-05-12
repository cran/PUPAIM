#'@title Weber-Van Vliet Isotherm Nonlinear Analysis
#'@name webervanvlietanalysis
#'@description It provides an excellent description of data patterns for a broad
#'range of systems. This model is suitable for batch rate and fixed-bed modelling
#'procedures as it gives a direct parameter evaluation.
#'@param Ce the numerical value for the equilibrium capacity
#'@param Qe the numerical value for the adsorbed capacity
#'@import nls2
#'@import Metrics
#'@import stats
#'@import ggplot2
#'@return the nonlinear regression and the parameters for Weber-Van-Vliet
#'Isotherm Analysis
#'@examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#'@examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#'@examples webervanvlietanalysis(Qe,Ce)
#'@author Keith T. Ostan
#'@author Chester C. Deocaris
#'@references Van Vliet, B.M., Weber Jr., Hozumi, H.. (1979) <doi:10.1016/0043-1354(80)90107-4> Modeling and
#'prediction of specific compound adsorption by activated carbon and synthetic
#'adsorbents. Water Research Vol.14, pp. 1719 to 1728.
#'@export
#'

# Building the Weber-Van Vliet isotherm nonlinear form
webervanvlietanalysis<- function(Qe,Ce) {

  x <- Qe
  y <- Ce
  data <- data.frame(Qe,Ce)

# Weber-Van Vliet isotherm nonlinear equation
  fit1 <- y ~ P* x^(R*(x^s)+t)

# Setting of starting values
  N <- nrow(na.omit(data))
  start1 <- data.frame(P = seq(0, 100, length.out = N),
                      R = seq(-1, 10, length.out = N),
                      s = seq(-1, 1, length.out = N),
                      t = seq(-1, 1, length.out = N))

# Fitting of the Weber-Van Vliet isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 1000, warnOnly = TRUE),
                 algorithm = "port")

  print("Weber Van-Vliet Isotherm Non-linear Analysis")
  print(summary(fit2))

  print("Aikake Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Weber-Van Vliet isotherm model

  error <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Predicted Values",
         "Relative Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- error(y)
  print(a)

# Graphical representation of the Weber-Van Vliet isotherm model

  ### Predicted parameter values
  parsweberV <- as.vector(coefficients(fit2))
  pars_P <- parsweberV[1L];
  pars_R <- parsweberV[2L];
  pars_s <- parsweberV[3L];
  pars_t <- parsweberV[4L];

  rhs <- function(x){(pars_P* x^(pars_R*(x^pars_s)+pars_t))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Qe",
         y = "Ce",
         title = "Weber-Van Vliet Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5)) + ggplot2::coord_flip()
}
