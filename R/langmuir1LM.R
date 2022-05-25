#' @title Langmuir Isotherm First Linear Form Analysis
#' @name langmuir1.LM
#' @description The Langmuir isotherm is described to be the most useful and
#' simplest isotherm for both chemical adsorption and physical adsorption. It
#' assumes that there is uniform adsorption energy onto the monolayer surface
#' and that there would be no interaction between the adsorbate and the surface.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the parameters for the Langmuir isotherm (first form), model error analysis,
#' and linear regression analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples langmuir1.LM(Ce,Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Langmuir, I. (1918) <doi:/10.1021/ja01269a066> The adsorption of
#' gases on plane surfaces of glass, mics and platinum. Journal of the American
#' Chemical Society, 1361-1403.
#' @references Chen, X. (2015) <doi:/10.3390/info6010014> Modeling of Experimental
#' Adsorption Isotherm Data. 14-22.
#' @export

# Building the Langmuir isotherm linear form
langmuir1.LM <- function(Ce,Qe) {
  x <- Ce
  y <- Ce/Qe
  data <- data.frame(x,y)

  # Fitting of the Langmuir isotherm linear form
  rhs <- function (x, Kl, Qmax) {
    (1/(Kl * Qmax))  + (x/Qmax)
  }

  fit1 <- lm(y ~ x)

  print("Langmuir Isotherm First Linear Form Analysis")
  print(summary(fit1))
  a <- (summary(fit1))

  ### Parameter values calculation
  qmax <- 1/a$coefficients[2]
  print("Qmax")
  print(qmax)

  kl <- ((1/a$coefficients[1])*(1/qmax))
  print("Kl")
  print(kl)
# -------------------------------------------

  print("Aikake Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

# Error analysis of the Langmuir isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae(y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data))
    SE <- SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Root Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Langmuir isotherm model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "Ce",
                  y = "Ce/Qe",
                  title = "Langmuir Isotherm Linear Model",
                  subtitle = "First Linear Form",
                  caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5))
}
