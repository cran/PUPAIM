#' @title Fowler-Guggenheim Isotherm Non-Linear Analysis
#' @name fowlerguggenheimanalysis
#' @description In Fowler-Guggenheim isotherm model, the lateral interaction of
#' the adsorbed molecules is taken into consideration. This is formulated on the
#'  basis that the heat adsorption process may vary positively or negatively
#'  with loading.
#' @param Ce is equal to Co which is the numeric value for the initial concentration
#' @param theta is the fractional surface coverage
#' @param Temp temperature
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Fowler-Guggenheim isotherm,
#'  and model error analysis
#' @examples theta <- c(0.19729, 0.34870, 0.61475, 0.74324, 0.88544, 0.89007, 0.91067, 0.91067, 0.96114)
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Temp <- 298
#' @examples fowlerguggenheimanalysis(theta,Ce,Temp)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Fowler, R. H. and Guggenheim, E. A. (1939) Statistical
#' Thermodynamics, Cambridge University Press, London, England.
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export
#'

# Building the Fowler-Guggenheim isotherm nonlinear form
fowlerguggenheimanalysis <- function(theta, Ce, Temp){

  x <- theta
  y <- Ce
  t <- Temp
  R <- 8.314
  data <- data.frame(x, y)


#  Fowler-Guggenheim nonlinear equation
  fit1 <- (y) ~ (1/KFG)*((x/(1-x))*exp((2*W*x)/(R*t)))

# Setting of starting values
  start1 <- data.frame(W = c(1,100), KFG = c(1,100))

# Fitting of the Fowler-Guggenheim isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1,  data=data,
                control = nls.control(maxiter = 100, warnOnly = TRUE),
                algorithm = "default")

  print("Fowler Guggenheim Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Fowler-Guggenheim isotherm model

  errors <- function(y){
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae (y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Square Error"= rmse,
         "Mean Absolute Error"= mae,
         "Mean Squared Error"= mse,
         "Relative Absolute Error"= rae,
         "Standard Error for the Regression S" = SE)
  }
  hga <- errors(y)
  print(hga)

  # Graphical representation of the Fowler-Guggenheim isotherm model

  ### Predicted parameter values
  parsFowlerG <- as.vector(coefficients(fit2))
  pars_w <- parsFowlerG[1L];
  pars_KFG <- parsFowlerG[2L];

  rhs <- function(x){(1/pars_KFG)*((x/(1-x))*exp((2*pars_w*x)/(R*t)))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = expression(paste(theta)),
         y = "Qe",
         title = "Fowler-Guggenheim Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

