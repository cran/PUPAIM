#' @title Dubinin-Radushkevich Isotherm  Linear Analysis
#' @name dubininraduskevich.LM
#' @description Dubinin-Radushkevich isotherm model is being utilized to define
#' adsorption energy mechanisms with Gaussian distribution onto heterogeneous surfaces.
#' Specifically, this model works well with an intermediate range of adsorbate
#' concentrations because it shows abnormal asymptotic behavior and is unable to
#' forecast Henry's Law at low pressure.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @param Temp temperature
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Dubinin-Radushkevich isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples Temp <- 298
#' @examples dubininradushkevich.LM (Ce,Qe,Temp)
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


# Building the Dubinin-Radushkevich isotherm linear form
dubininradushkevich.LM <- function(Ce, Qe, Temp){

  t <- Temp
  R <- 8.314
  epsilon <- R*t*log(1+(1/Ce))
  x <- epsilon^2
  y <- log(Qe)
  data <- data.frame(x,y)

  # Fitting of the Dubinin-Radushkevich isotherm linear form
  rhs <- function(x,qs,K) {
    log(qs) - K*(epsilon^2)
  }

  fit1 <- lm(y~x)

  print("Dubinin-Radushkevich Isotherm Linear Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- (summary(fit1))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  qs <- exp(a)
  print("qs")
  print(qs)

  K <- -b
  print("K")
  print(K)
  # ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

  # Error analysis of the Dubinin-Radushkevich isotherm linear model

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

# Graphical representation of the Dubinin-Radushkevich isotherm linear model

#### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste(epsilon^"2")),
         y = "ln(Qe)",
         title = "Dubinin-Radushkevich Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
