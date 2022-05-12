#' @title Elovich Isotherm Linear Analysis
#' @name elovich.LM
#' @description Elovich isotherm model is based on kinetic principle which
#' assumes that the adsorption sites would exponentially increase with chemical
#' reactions responsible for adsorption. It is suited for describing the behavior
#' of adsorption concurring with the nature of chemisorption.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Elovich isotherm, and model
#' error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples elovich.LM(Ce,Qe)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Zeldowitsch, J. (1934). "über Den Mechanismus der Katalytischen
#' Oxidation Von CO a MnO2," URSS, Acta Physiochim, Vol. 1, No. 2, 1934, pp. 364-449.
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export

# Building the Elovich isotherm linear form
elovich.LM <- function(Ce,Qe){

  x <- Qe
  y <- log(Qe/Ce)
  data <- data.frame(x, y)

# Fitting of the Elovich isotherm linear form
   rhs <- function(x, KE, Qmax){
    log(KE*Qmax)-(x/Qmax)
  }

  fit1 <- lm(y~x)

  print("Elovich Analysis")
  print(summary(fit1))

  ### y = a + bx
  c <- summary(fit1)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  Qmax <- -1/b
  print("Qmax")
  print(Qmax)

  KE <- exp(a)/Qmax
  print("KE")
  print(KE)
# ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

# Error analysis of the Elovich iostherm model
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
         "Relative Absolute Error" = rae,
         "Standard Error of the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Elovich isotherm linear model

#### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "Qe",
         y = "ln(Qe/Ce)",
         title = "Elovich Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
