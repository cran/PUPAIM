#' @title Hill Isotherm Linear Analysis
#' @name hill.LM
#' @description Hill isotherm model shows the connection of different species
#' being adsorbed on to the homogeneous surfaces. This isotherm model supposes
#' that adsorption is a cooperative phenomenon which means the adsorbates having
#' the capability to bind at one specific site on the adsorbent affecting other
#' binding sites on the same adsorbent
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for the Hill isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples hill.LM(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Hill, T. L. (1946). <doi:10.1063/1.1724129> "Statistical mechanics of
#' multimolecular adsorption II. Localized and mobile adsorption and absorption,"
#' The Journal of Chemical Physics, vol. 14, no. 7, pp. 441-453.
#' @export

# Building the Hill isotherm linear form
hill.LM<- function(Ce, Qe){

  x1 <- Ce
  y1 <- Qe
  data <- data.frame(x1, y1)

# Hill isotherm nonlinear equation
  fit1 <- y1 ~ ((x1^nh)/(Kd+x1^nh))

# Setting of starting values
  start1 <- list(nh = 1, Kd = 1)

# Fitting of the Hill isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
               control = nls.control(maxiter = 100, warnOnly = TRUE),
               algorithm = "plinear")

  param <- summary(fit2)
  qh <- param$coefficients[3]

# Establishing Hill isotherm linear form
  x <- log10(Ce)
  y <- log10(Qe/(qh-Qe))

  rhs1 <- function (x, Kd, nh) {
    nh*log10(x)-log10(Kd)
  }

# Hill isotherm linear fitting
  fit3 <- lm(y~x)

  print("Hill Isotherm Linear Analysis")
  print(summary(fit3))

  ### y = a + bx
  c <- summary(fit3)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  nh <- b
  print("nh")
  print(nh)

  Kd <- -a
  print("Kd")
  print(Kd)

  print("qh")
  print(qh)
# -------------------------------------
  print("Aikake Information Criterion")
  print(AIC(fit3))

  print("Bayesian Information Criterion")
  print(BIC(fit3))

# Error analysis of the Hill isotherm linear model
  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit3))
    mae <- Metrics::mae(y, predict(fit3))
    mse <- Metrics::mse(y, predict(fit3))
    rae <- Metrics::rae(y, predict(fit3))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit3))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Hill isotherm linear model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "log(Ce)",
         y = "log(Qe/(Qsh-Qe))",
         title = "Hill Isotherm Linear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
