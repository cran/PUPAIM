#' @title Kiselev Isotherm Linear Analysis
#' @name kiselev.LM
#' @description It is also known as localized monomolecular layer model and
#' is only valid for surface coverage theta > 0.68.
#' @param Ce the numerical value for equilibrium capacity
#' @param theta is the fractional surface coverage
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for the Kiselev isotherm, and
#' model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607,
#' 0.80435, 1.10327, 1.58223)
#' @examples theta  <- c(0.1972984, 0.3487013, 0.6147560, 0.7432401, 0.8854408,
#' 0.8900708, 0.9106746, 0.9106746, 0.9611422)
#' @examples kiselev.LM(theta, Ce)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Kiselev, A. V. (1958). "Vapor adsorption in the formation of
#' adsorbate molecule complexes on the surface," Kolloid Zhur, vol. 20, pp. 338-348.
#' @export

# Building the Kiselev isotherm linear form
kiselev.LM <- function(theta,Ce){

  x <- 1/theta
  y <- 1/(Ce*(1-theta))
  data <- data.frame(x, y)

# Fitting of the Kiselev isotherm linear form
  rhs <- function(x, Ki, Kn) {
  (Ki * Kn) + (Ki/theta)
  }

  fit1 <- lm(y~x)

  print("Kiselev Isotherm Linear Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- summary(fit1)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  Ki <- b
  print("Ki")
  print(Ki)

  Kn <- a/b
  print("Kn")
  print(Kn)
  # ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

# Error analysis of the Kiselev isotherm model
  errors <- function(y){
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae (y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y)<- rownames(theta) <- colnames(theta)
    list("Relative Mean Square Error"= rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error"= mse,
         "Relative Absolute Error"= rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Kiselev isotherm model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste("1/", theta)),
         y = expression(paste("1/Ce(1-", theta,")")),
         title = "Kiselev Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}


