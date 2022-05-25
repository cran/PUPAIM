#' @title Temkin Isotherm Linear Analysis
#' @name temkin.LM
#' @description Temkin isotherm  is a monolayer adsorption isotherm model which
#' takes into account the effects that the indirect interaction amongst adsorbate
#' molecules could have on the adsorption process.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @param Temp temperature
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Temkin isotherm, and model
#' error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples Temp <- 298.15
#' @examples temkin.LM(Ce,Qe,Temp)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Temkin, M.J., and Pyzhev, V. (1940). Kinetics of ammonia synthesis
#' on promoted iron catalyst. Acta Phys. Chim. USSR 12, 327-356.
#' @references Foo, K. Y., and Hameed, B. H. (2009, September 13). <doi:10.1016/j.cej.2009.09.013>
#' Insights into the modeling of adsorption isotherm systems. Chemical Engineering Journal.
#' @export
#'

# Building the Temkin isotherm linear form
temkin.LM <- function(Ce, Qe, Temp){

  x <- log(Ce)
  y <- Qe
  t <- Temp
  R <- 8.314
  data <- data.frame(x,y)

# Fitting of the Temkin isotherm linear form
rhs <- function(x,aT,bT) {
     ((R*t)/bT)*log(aT) + ((R*t)/bT)* log(x)
    }

  fit1 <- lm(y~x)

  print("Temkin Isotherm Linear Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- summary(fit1)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  bT <- (b/(R*t))^-1
  print("bT")
  print(bT)

  aT <- exp(a*(bT/(R*t)))
  print("aT")
  print(aT)
# ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

# Error analysis of the Temkin isotherm linear model
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

# Graphical representation of the Temkin isotherm linear model

#### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "ln(Ce)",
         y = "Qe",
         title = "Temkin Isotherm Linear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
