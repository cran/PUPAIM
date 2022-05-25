#' @title Sips Isotherm Linear Analysis
#' @name sips.LM
#' @description It is the most applicable to use in the monolayer adsorption
#' isotherm model amongst the three-parameter isotherm models and is also valid
#' for the prediction of heterogeneous adsorption systems as well as localized
#' adsorption with no interactions occurring between adsorbates.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Sips isotherm, and model error
#' analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples sips.LM(Ce,Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Sips, R. (1948) <doi:10.1063/1.1746922> On the structure of a catalyst
#' surface. The Journal of Chemical Physics, 16(5), 490-495.
#' @export
#'

# Building the Sips isotherm linear form
sips.LM <- function(Ce, Qe) {

  x1 <- Ce
  y1 <- Qe
  data<- data.frame(x1, y1)

# Obtaining the model exponent
### Sips isotherm nonlinear equation
  fit1 <- y1 ~ (Ks*(x1^Ns))/(1 + (As*(x1^Ns)))

### Setting of starting values
  start1 <- list(As = 1, Ks = 1, Ns = 1)

### Fitting of Sips isotherm via nls2
  fit2 <- nls2::nls2(fit1, start = start1, data=data,
               control = nls.control(maxiter = 50, warnOnly = TRUE),
               algorithm = "port")

  param <- summary(fit2)
  expModel <- param$coefficients[3]

# Establishing Sips isotherm linear form
  x <- 1/Ce^expModel
  y <- 1/Qe

  rhs <- function (x, As, Ks) {
    (As/Ks) + (1/Ks)*x
    }

# Sips isotherm linear fitting
  fit3 <- lm(y~x)

  print("Sips Isotherm Linear Analysis")
  print(summary(fit3))

### y = a + bx
  c <- (summary(fit3))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

### Parameter values calculation
  Ks <- b^-1
  print("Ks")
  print(Ks)

  As <- a*Ks
  print("As")
  print(As)

  Ns <- expModel
  print("Ns")
  print(Ns)
# -------------------------------------
  AIC <- AIC(fit3)
  print("Aikake Information Criterion")
  print(AIC)

  BIC <- BIC(fit3)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Sips isotherm linear model
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

# Graphical representation of the Sips isotherm linear model

#### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste("Ce"^"Ns")),
         y = "1/Qe",
         title = "Sips Isotherm Linear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}

