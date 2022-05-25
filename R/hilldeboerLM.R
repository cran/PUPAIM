#' @title Hill-Deboer Isotherm Linear Analysis
#' @name hilldeboer.LM
#' @description Hill-Deboer isotherm model describes as a case where there is
#' mobile adsorption as well as lateral interaction among molecules.The increased
#' or decreased affinity depends on the kind of force among the adsorption molecules.
#' If there is an attraction between adsorbed molecules, there is an increase in
#' affinity. On the other hand, decreased affinity happens when there is repulsion
#' among the adsorbed molecules.
#' @param Ce the numerical value for the equilibrium capacity
#' @param theta is the fractional surface coverage
#' @param Temp temperature
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for the Hill-Deboer isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607,
#' 0.80435, 1.10327, 1.58223)
#' @examples theta  <- c(0.1972984, 0.3487013, 0.6147560, 0.7432401, 0.8854408,
#' 0.8900708, 0.9106746, 0.9106746, 0.9611422)
#' @examples Temp <- 298.15
#' @examples hilldeboer.LM(Ce,theta, Temp)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references De Boer, J. H. (1953). The Dynamical Character of adsorption,
#' Oxford University Press, Oxford, England.
#' @references Foo, K. Y., and Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export

# Building the Hill-Deboer isotherm linear forms
hilldeboer.LM <- function(Ce,theta, Temp){

  x <- theta
  y <- log((Ce*(1-theta))/theta) - (theta/(1-theta))
  R <- 8.314
  t <- Temp
  data <- data.frame(x, y)

  # Fitting of the Hill-Deboer isotherm linear form
  rhs <- function(x, K1, K2) {
    -log(K1) - (K2 - theta)/(R*t)
  }

  fit1 <- lm(y~x)

  print("Hill-Deboer Isotherm Linear Analysis")
  print(summary(fit1))

  ### y = a+bx
  c <- (summary(fit1))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  K1 <- -exp(a)
  print("K1")
  print(K1)

  K2 <- -b*(R*t)
  print("K2")
  print(K2)

  # ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

  # Error analysis of the Hill-Deboer isotherm model
  errors <- function(y){
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae (y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y)<- rownames(y) <- colnames(y)
    list("Relative Mean Square Error"= rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error"= mse,
         "Relative Absolute Error"= rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  # Graphical representation of the Hill-Deboer isotherm model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x , method = "lm", se = F, color = "#D35400") +
    ggplot2::labs(x = expression(paste(theta)),
         y = expression(paste("log((Ce(1-",theta,"))/", theta,")-(",theta,"/(1-",theta,")")),
         title = "Hill-de Boer Isotherm Linear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
