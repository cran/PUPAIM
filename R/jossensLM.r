#' @title Jossens Isotherm Linear Analysis
#' @name jossens.LM
#' @description The Jossens isotherm model predicts a simple equation based on
#' the energy distribution of adsorbate-adsorbent interactions at adsorption
#' sites. This model assumes that the adsorbent has heterogeneous surface with
#' respect to the interactions it has with the adsorbate.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for the Jossens isotherm, and model error
#' analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples jossens.LM(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Jossens, L., Prausnitz, J. M., Fritz, W., Schlunder, E. U.,
#' & Myers, A. L. (1978) <doi:10.1016/0009-2509(78)85015-5> Thermodynamics of
#' multi-solute adsorption from dilute aqueous solutions.
#' Chemical Engineering Science, 33(8), 1097-1106.
#' @export
#'

# Building the Jossens isotherm linear form
jossens.LM <- function(Ce,Qe) {

  x1 <- Qe
  y1 <- Ce
  data <- data.frame(x1, y1)

  # Obtaining the model exponent
  ### Jossens isotherm nonlinear equation
  fit1 <- (y1 ~ (x1/H)*(exp(J*(x1^Nj))))

  ### Setting of starting values
  start1 <- data.frame(H = c(1, 500), J = c(1, 100), Nj = c(0, 1))

  ### Fitting of the Jossens isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1,data=data,
               control = nls.control(maxiter = 100 , warnOnly = TRUE),
               algorithm = "port")

  param <- summary(fit2)
  expModel <- param$coefficients[3]

  # Establishing Jossens isotherm linear form
  x <- Ce^expModel
  y <- log(Ce/Qe)

  rhs <- function (x, H, J) {
    -log(H) + J*x
  }

  # Jossens isotherm linear fitting
  fit3 <- lm(y~x)

  print("Jossens Isotherm Linear Analysis")
  print(summary(fit3))

  ### y = a + bx
  c <- summary(fit3)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  J <- b
  print("J")
  print(J)

  H <- -exp(a)
  print("H")
  print(H)

  Nj <- expModel
  print("Nj")
  print(Nj)
  # -------------------------------------
  print("Aikake Information Criterion")
  print(AIC(fit3))

  print("Bayesian Information Criterion")
  print(BIC(fit3))

  # Error analysis of the Jossens isotherm linear model
  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit3))
    mae <- Metrics::mae(y, predict(fit3))
    mse <- Metrics::mse(y, predict(fit3))
    rae <- Metrics::rae(y, predict(fit3))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a<- errors(y)
  print(a)

  # Graphical representation of the Jossens isotherm linear model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste("Qe"^"Nj")),
         y ="ln(Ce/Qe)",
         title = "Jossens Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}

