#' @title Halsey Isotherm Linear Analysis
#' @name halsey.LM
#' @description A multilayer adsorption isotherm model which is suited for
#' adsorption of adsorbate ions at a distance that is relatively large from the
#' surface.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for the Halsey isotherm, and
#' model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples halsey.LM(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Halsey, G., & Taylor, H. S. (1947) <doi:10.1063/1.1746618>
#' The adsorption of hydrogen on tungsten powders. The Journal of
#' Chemical Physics, 15(9), 624-630.
#' @export
#'

# Building the Halsey isotherm linear form
halsey.LM <- function(Ce, Qe){

  x <- log(Ce)
  y <- log(Qe)
  data <- data.frame(x, y)

# Halsey isotherm linear equation
  rhs <- function(x, Kh, nh) {
    (((1/nh) * log(Kh)) - ((1/nh) * log(x)))
  }

# Fitting of the Halsey isotherm
  fit1 <- lm(y~x)

  print("Halsey Isotherm Linear Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- (summary(fit1))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  nh <- (-b)^-1
  print("nh")
  print(nh)

  Kh <- exp(a*nh)
  print("Kh")
  print(Kh)
# ---------------------------------
  AIC <- AIC(fit1)
  print("Akaike Information Criterion")
  print(AIC)

  BIC <- BIC(fit1)
  print("Bayesian Information Criterion")
  print(BIC)

# Error Analysis of the Halsey isotherm model
  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae(y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
    }
  a <- errors(y)
  print(a)

# Graphical representation of the Halsey isotherm model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "ln(Ce)",
         y = "ln(Qe)",
         title = "Halsey Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
