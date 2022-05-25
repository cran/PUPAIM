#' @title Henry Isotherm Linear Analysis
#' @name henryanalysis
#' @description It describes the appropriate fit to the adsorption of adsorbate
#' at relatively low concentrations such that all adsorbate molecules are
#' secluded from their nearest neighbours.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the linear regression, parameters for the Henry isotherm, and model
#' error analysis
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples henryanalysis(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Deocaris, C., and Osio, L. (2020). Fitting Henry's
#' Adsorption Isotherm model in R using PUPAIM package.
#' @export

# Building the Henry isotherm linear model
henryanalysis <- function (Ce, Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x,y)

# Henry isotherm linear equation
  rhs <- function(x, K){
  (K*x)
    }

# Fitting of Henry isotherm
  fit230 <- lm(y ~ x)

  print("Henry Isotherm Analysis")
  print(summary(fit230))

  a <- (summary(fit230))
  b <- a$coefficients[2]

  ### Parameter values calculation
  print("K")
  print(b)

# ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit230))

  print("Bayesian Information Criterion")
  print(BIC(fit230))

# Error analysis of the Henry isotherm linear model
  errors <- function(y){
    rmse <- Metrics::rmse(y, predict(fit230))
    mae <- Metrics::mae (y, predict(fit230))
    mse <- Metrics::mse(y, predict(fit230))
    rae <- Metrics::rae(y, predict(fit230))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit230))^2)/(N-2))
    colnames(y)<- rownames(y) <- colnames(y)
    list("Relative Mean Square Error"= rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error"= mse,
         "Relative Absolute Error"= rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Henry isotherm model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Henry Isotherm Linear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
