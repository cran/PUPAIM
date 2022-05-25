#' @title Harkins-Jura Isotherm Non-Linear Analysis
#' @name harkinsjuraanalysis
#' @description A model that assumes the possibility of multilayer adsorption
#' on the surface of absorbents having heterogenous pore distribution
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Harkins-Jura isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples harkinsjuraanalysis(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Harkins, W. D., and Jura, G. (1944) <doi:10.1021/ja01236a048> Surfaces of solids. XIII.
#' A vapor adsorption method for the determination of the area of a solid without
#' the assumption of a molecular area, and the areas occupied by nitrogen and other
#' molecules on the surface of a solid. Journal of the American Chemical Society,
#' 66(8), 1366-1373.
#' @export

# Building the Harkins-Jura isotherm nonlinear form
harkinsjuraanalysis <- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x,y)

# Harkins-Jura isotherm nonlinear equation
  fit1 <- y ~ (A/(b-log(x)))^1/2

# Setting of starting values
  start1 <- data.frame(A = c(1, 100), b = c(1, 100))

# Fitting of the Harkins-Jura isotherm via nls2
  fit2 <- nls2::nls2(fit1, start = start1,  data=data,
                control = nls.control(maxiter = 100 , warnOnly = TRUE),
                algorithm = "port")

  print("Harkins-Jura Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Harkins-Jura isotherm model
 errors <- function(y) {
   rmse <- Metrics::rmse(y, predict(fit2))
   mae <- Metrics::mae(y, predict(fit2))
   mse <- Metrics::mse(y, predict(fit2))
   rae <- Metrics::rae(y, predict(fit2))
   N <- nrow(na.omit(data))
   SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
   colnames(y) <- rownames(y) <- colnames(y)
   list("Root Mean Squared Error" = rmse,
        "Mean Absolute Error" = mae,
       "Mean Squared Error" = mse,
       "Relative Absolute Error" = rae,
       "Standard Error for the Regression S" = SE)
  }
  a <- errors (y)
  print(a)

  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))

# Graphical representation of the Harkins-Jura isotherm model

 ### Predicted parameter values
 parsharkj <- as.vector(coefficients(fit2))
 pars_A <- parsharkj[1L];
 pars_b <- parsharkj[2L]

 rhs <- function(x) {(pars_A/(pars_b-log(x)))^1/2}

 #### Plot details
 ggplot2::theme_set(ggplot2::theme_bw(10))
 ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
   ggplot2::geom_function(color = "#D35400", fun = rhs ) +
   ggplot2::labs(x = "Ce",
        y = "Qe",
        title = "Harkins-Jura Isotherm Nonlinear Model",
        caption = "PUPAIM") +
   ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
