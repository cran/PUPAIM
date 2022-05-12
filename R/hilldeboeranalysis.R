#' @title Hill-Deboer Isotherm Non-Linear Analysis
#' @name hilldeboeranalysis
#' @description Hill-Deboer isotherm model describes as a case where there is
#' mobile adsorption as well as lateral interaction among molecules.The increased
#' or decreased affinity depends on the kind of force among the adsorption molecules.
#' If there is an attraction between adsorbed molecules, there is an increase in
#' affinity. On the other hand, decreased affinity happens when there is repulsion
#' among the adsorbed molecules.
#' @param theta is the fractional surface coverage
#' @param Ce the numerical value for the equilibrium capacity
#' @param Temp the temperature of the adsorption experimentation in Kelvin
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Hill-Deboer isotherm,
#'and model error analysis
#' @examples theta <- c(0.19729, 0.34870, 0.61475, 0.74324, 0.88544, 0.89007, 0.91067, 0.91067, 0.96114)
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Temp <- 298
#' @examples hilldeboeranalysis(theta, Ce, Temp)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references De Boer, J. H. (1953). The Dynamical Character of adsorption,
#' Oxford University Press, Oxford, England.
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export

# Building the Hill-De Boer isotherm nonlinear forms
hilldeboeranalysis <- function(theta, Ce, Temp){

  x <- theta
  y <- Ce
  t <- Temp
  R <- 8.314
  data <- data.frame(x, y)

# Hill-De Boer isotherm nonlinear equation
  fit1 <- y ~ (x/(K1*(1-x)))*exp((x/(1-x))-((K2*x)/R*t))

# Setting of starting values
  N <- nrow(na.omit(data))
  start1 <- data.frame(K1 = seq(1, 100, length.out = N),
                       K2 = seq(1, 100, length.out = N))

# Fitting of the Hill-Deboer isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 100, warnOnly = TRUE),
                 algorithm = "port")

  print("Hill-Deboer Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Hill-Deboer isotherm model

errors <- function(y){
  rmse <- Metrics::rmse(y, predict(fit2))
  mae <- Metrics::mae (y, predict(fit2))
  mse <- Metrics::mse(y, predict(fit2))
  rae <- Metrics::rae(y, predict(fit2))
  SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
  colnames(y) <- rownames(y) <- colnames(y)
  list('Root Mean Square Error'= rmse,
       'Mean Absolute Error'= mae,
      'Mean Squared Error'= mse,
      'Relative Absolute Error'= rae,
      'Standard Error for the Regression S' = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Hill-Deboer isotherm model

  ### Predicted parameter values
  parshilldeb <- as.vector(coefficients(fit2))
  pars_K1 <- parshilldeb[1L];
  pars_K2 <- parshilldeb[2L];

  rhs <- function(x){((x/(pars_K1*(1-x)))*exp((x/(1-x))-((pars_K2 *x)/R*t)))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = expression(paste(theta)),
         y = "Ce",
         title = "Hill-de Boer Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
