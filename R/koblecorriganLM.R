#' @title Koble-Carrigan Isotherm Linear Analysis
#' @name koblecarrigan.LM
#' @description It is three-parameter isotherm model equation that incorporates
#' both Freundlich and Langmuir isotherms for representing equilibrium adsorption
#' data. Koble-Corrigan isotherm model appeared to have advantages over both the
#' Langmuir and Freundlich equations in that it expresses adsorption data over
#' very wide ranges of pressures and temperatures.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Koble-Carrigan isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples koblecarrigan.LM(Ce, Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Corrigan, T. E., & Koble, R. A.(1952) <doi:10.1021/ie50506a049>
#' Adsorption isotherms for pure hydrocarbons Ind. Eng. Chem. 44 383-387.
#' @export

# Building the Sips isotherm linear form
koblecarrigan.LM <- function(Ce,Qe) {

  x1 <- Ce
  y1 <- Qe
  data <- data.frame(x1, y1)

# Koble-Corrigan isotherm nonlinear equation
  fit1 <- y1 ~ (Akc*(x1^Nkc))/(1 + Bkc*(x1^Nkc))

# Setting of starting values
  start1 <- list(Akc = 1, Bkc = 1, Nkc = 1)

# Fitting of the Koble-Corrigan isotherm via nls2
  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter =50 , warnOnly = TRUE),
                 algorithm = "port")

  param <- summary(fit2)
  expModel <- param$coefficients[3]

# Establishing Koble-Corrigan isotherm linear form
  x <- 1/(Ce^expModel)
  y <- 1/Qe

  rhs <- function (x, Bkc, Akc) {
    (Bkc/Akc) + (1/Akc*x)
  }

# Koble-Corrigan  isotherm linear fitting
  fit3 <- lm(y~x)

  print("Koble-Corrigan  Isotherm Linear Analysis")
  print(summary(fit3))

### y = a + bx
  c <- (summary(fit3))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

### Parameter values calculation
  Akc <- b^-1
  print("Akc")
  print(Akc)

  Bkc <- a*Akc
  print("Bkc")
  print(Bkc)

  Nkc <- expModel
  print("Nkc")
  print(Nkc)
# -------------------------------------
  AIC <- AIC(fit3)
  print("Aikake Information Criterion")
  print(AIC)

  BIC <- BIC(fit3)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Koble-Corrigan isotherm linear model
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

# Graphical representation of the Koble-Corrigan isotherm linear model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste("1"/"Ce"^"Nkc")),
         y = "1/Qe",
         title = "Koble-Corrigan Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
