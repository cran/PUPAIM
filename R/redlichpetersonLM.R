#' @title Redlich-Peterson Isotherm Linear Analysis
#' @name redlichpeterson.LM
#' @description Redlich-Peterson isotherm model has an exponential function
#' which can be found in the denominator and in the numerator, it has a linear
#' dependence on the concentration denoting the adsorption equilibrium depending
#' on a wide range of concentration
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Redlich-Peterson isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples redlichpeterson.LM(Ce,Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Peterson, D. L. and Redlich, O.(1959) <doi:10.1021/j150576a611>
#' A useful adsorption isotherm. J PhysChem US;63(6):1024. Research, vol. 6, no. 1,
#' pp. 265-276, 2012.
#' @export
#'

# Building the Redlich-Peterson isotherm linear form
redlichpeterson.LM <- function(Ce,Qe) {

  x1 <- Ce
  y1 <- Qe
  data <- data.frame(x1, y1)

# Obtaining the model exponent
### Redlich-Peterson isotherm nonlinear equation
  fit1 <- y1 ~ (Krep*x1)/(1+(Arep*(x1^Nrep)))

### Setting of starting values
  start1 <- list(Krep = 1, Arep = 1, Nrep = 1)

### Fitting of the Redlich-Peterson isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 100 , warnOnly = TRUE),
                 algorithm = "port")

  param <- summary(fit2)
  expModel <- param$coefficients[3]

# Establishing Redlich-Peterson isotherm linear form
  x <- Ce^expModel
  y <- Ce/Qe

  rhs <- function (x, Qmax, Krp) {
    1/(Krep) + (Arep/Krep)*x
  }

# Redlich-Peterson isotherm linear fitting
  fit3 <- lm(y~x)

  print("Redlich-Peterson Isotherm Linear Analysis")
  print(summary(fit3))

### y = a + bx
  c <- (summary(fit3))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

### Parameter values calculation
  Krep <- a^-1
  print("Krep")
  print(Krep)

  Arep <- b*Krep
  print("Arep")
  print(Arep)

  Nrep <- expModel
  print("Nrep")
  print(Nrep)
# -------------------------------------
  AIC <- AIC(fit3)
  print("Aikake Information Criterion")
  print(AIC)

  BIC <- BIC(fit3)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Redlich-Peterson isotherm linear model
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
  a<- errors(y)
  print(a)

# Graphical representation of the Redlich-Peterson isotherm linear model

#### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste("Ce"^"g")),
         y = "Ce/Qe",
         title = "Redlich-Peterson Isotherm Linear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
