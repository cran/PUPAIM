#' @title Fowler-Guggenheim Isotherm Linear Analysis
#' @name fowlerguggenheim.LM
#' @description In Fowler-Guggenheim isotherm model, the lateral interaction of
#' the adsorbed molecules is taken into consideration. This is formulated on the
#' basis that the heat adsorption process may vary positively or negatively with
#' loading.
#' @param Ce is equal to the numerical value for the equilibrium capacity
#' @param theta is fractional surface coverage
#' @param Temp temperature
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Fowler-Guggenheim isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607,
#' 0.80435, 1.10327, 1.58223)
#' @examples theta  <- c(0.1972984, 0.3487013, 0.6147560, 0.7432401, 0.8854408,
#' 0.8900708, 0.9106746, 0.9106746, 0.9611422)
#' @examples Temp <- 298.15
#' @examples fowlerguggenheim.LM(theta, Ce, Temp)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Fowler, R. H. and Guggenheim, E. A. (1939) Statistical Thermodynamics,
#' Cambridge University Press, London, England.
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export
#'

# Building the Fowler-Guggenheim isotherm linear form
fowlerguggenheim.LM <- function(theta, Ce, Temp){

  x <- Ce
  y <- log((Ce*(1-theta))/theta)
  data <- data.frame(x, y)
  t <- Temp
  R <- 8.314

  # Fitting of the Fowler-Guggenheim isotherm linear form
  rhs <- function(theta, KFG, W){
    -log(KFG) + ((2*W*theta)/(R*t))
  }

  fit1<- lm(y~x)

  print("Fowler-Guggenheim Isotherm Linear Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- summary(fit1)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  KFG <- -exp(a)
  print("KFG")
  print(KFG)

  W <- ((R*t)/2)*b
  print("W")
  print(W)
# ---------------------------------
  AIC <- AIC(fit1)
  print("Akaike Information Criterion")
  print(AIC)

  BIC <- BIC(fit1)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Fowler-Guggenheim isotherm model
  errors <- function(y){
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae (y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Square Error"= rmse,
         "Mean Absolute Error"= mae,
         "Mean Squared Error"= mse,
         "Relative Absolute Error"= rae,
         "Standard Error for the Regression S" = SE)
  }
  z <- errors(y)
  print(z)

  # Graphical representation of the Fowler-Guggenheim isotherm linear model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste(theta)),
         y = expression(paste("ln(Ce(1-",theta,")/",theta,")")),
         title = "Fowler-Guggenheim Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
