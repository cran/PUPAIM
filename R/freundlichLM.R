#' @title Freundlich Isotherm Linear Analysis
#' @name freundlich.LM
#' @description This isotherm model is an empirical model applicable to diluted
#' solutions adsorption processes. Furthermore, this model gives an equation which
#' defines the surface heterogeneity and the exponential distribution of active sites.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Freundlich isotherm, and model
#' error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples freundlich.LM (Ce,Qe)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Freundlich, H. 1907. Ueber die adsorption in loesungen. Z.
#' Phys. Chem.57:385-470
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export
#'

# Building the Freundlich isotherm linear form
freundlich.LM <- function(Ce, Qe){

  x <- log10(Ce)
  y <- log10(Qe)
  data <- data.frame(x, y)

  # Fitting of the Freundlich isotherm linear form
  rhs <- function(x, Kf, n){
    log(Qe) ~ logKF + (1/n)*log(Ce)
  }

  fit1 <- lm(y~x)

  print("Freundlich Isotherm Linear Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- summary(fit1)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  Kf <- 10^(a)
  print("Kf")
  print(Kf)

  n <- 1/b
  print("n")
  print(n)
# ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

# Error analysis of the Freundlich isotherm model
  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae(y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Root Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Root Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  # Graphical representation of the Freundlich isotherm linear model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "log(Ce)",
         y = "log(Qe)",
         title = "Freundlich Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
