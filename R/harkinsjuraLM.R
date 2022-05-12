#' @title HarkinsJura Isotherm Linear Analysis
#' @name harkinsjura.LM
#' @description A model that assumes the possibility of multilayer adsorption
#' on the surface of absorbents having heterogenous pore distribution.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for the HarkinsJura isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples harkinsjura.LM(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Harkins, W. D., & Jura, G. (1944) <doi:10.1021/ja01236a048>
#' Surfaces of solids. XIII. A vapor adsorption method for the determination of the
#' area of a solid withoutthe assumption of a molecular area, and the areas occupied
#' by nitrogen and other molecules on the surface of a solid. Journal of the American
#' Chemical Society, 66(8), 1366-1373.
#' @export
#'

# Building the Harkins-Jura isotherm linear form
harkinsjura.LM<- function(Ce, Qe){

  x <- log10(Ce)
  y <- (1/(Qe)^2)
  data <- data.frame(x, y)

# Harkins-Jura isotherm linear equation
  rhs <- function(x, A, B) {
    B/A-(1/A)*log(x)
    }

# Fitting of the Harkins-Jura isotherm
  fit1 <- lm(y~x)

  print("Harkins-Jura Isotherm Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- summary(fit1)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  A <- (-b)^-1
  print("A")
  print(A)

  B <- a*A
  print("B")
  print(B)
# ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

# Error analysis of the Harkins-Jura isotherm model
errors <- function(y) {
  rmse <- Metrics::rmse(y, predict(fit1))
  mae <- Metrics::mae(y, predict(fit1))
  mse <- Metrics::mse(y, predict(fit1))
  rae <- Metrics::rae(y, predict(fit1))
  N <- nrow(na.omit(data))
  SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
  colnames(y) <- rownames(y) <- colnames(y)
  list("Relative Mean squared Error" = rmse,
       "Mean Absolute Error" = mae,
       "Mean Squared Error" = mse,
       "Relative Absolute Error" = rae,
       "Standard Error for the Regression S" = SE)
    }
  a <- errors(y)
  print(a)

# Graphical representation of the Harkins-Jura isotherm model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "log(Ce)",
         y = expression(paste("1/Qe"^"2")),
         title = "Harkins-Jura Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}


