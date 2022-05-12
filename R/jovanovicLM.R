#' @title Jovanovic Isotherm Linear Analysis
#' @name jovanovic.LM
#' @description The Jovanovic isotherm model was built upon the assumptions
#' based on the Langmuir isotherm model with few possible inclusions of mechanical
#' contact among the desorbing and adsorbing molecules. The adjustment of the
#' adsorption surface from this model made the equation less effective in the
#' physical adsorption but can be applied to adsorption with both mobile and
#' localized monolayer without lateral interaction. Moreover, the equation of
#' the Jovanovic isotherm model is able to reach the limit of saturation when
#' there is high concentration, while it reduces to Henry's Law at low concentration.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for the Jovanovic isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples jovanovic.LM(Ce,Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Jovanović, D.S. (1969) <doi:10.1007/BF01542531> Physical adsorption of gases.
#' Kolloid-Z.u.Z.Polymere 235, 1214–1225.
#' @export

# Building the Jovanovic isotherm linear form
jovanovic.LM <- function(Ce, Qe){

  x <- Ce
  y <- log(Qe)
  data <- data.frame(x,y)

  # Fitting of the Jovanovic isotherm linear form
  rhs <- function(x,Qmax) {
    y ~ log(Qmax) - Kj*Ce
  }

  fit1 <- lm(y~x)

  print("Jovanovic Isotherm Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- (summary(fit1))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  Qmax <- exp(a)
  print("Qmax")
  print(Qmax)

  Kj <- -b
  print("Kj")
  print(Kj)
  # ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

  # Error analysis of the Jovanovic isotherm linear model
  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae(y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean squared error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative absolute error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  # Graphical representation of the Jovanovic isotherm linear model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "Ce",
         y = "ln(Qe)",
         title = "Jovanovic Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
