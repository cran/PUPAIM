#' @title Halsey Isotherm Non-Linear Analysis
#' @name halseyanalysis
#' @description A multilayer adsorption isotherm model which is suited for
#' adsorption of adsorbate ions at a distance that is relatively large from the
#' surface.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for the Halsey isotherm, and
#' model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples halseyanalysis(Ce, Qe)
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Halsey, G., and Taylor, H. S. (1947) <doi:10.1063/1.1746618> The adsorption of
#' hydrogen on tungsten powders. The Journal of Chemical Physics, 15(9), 624-630.
#' @export

# Building the Halsey isotherm nonlinear form
halseyanalysis <- function(Ce,Qe) {

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# Halsey isotherm nonlinear equation
  fit1 <- y ~ exp((log(Kh)-log(x))/nh)

# Setting of starting values
  start1 <- list(nh = -1, Kh = 1)

# Fitting of the Halsey isotherm via nls2

  fit2 <- nls2::nls2(fit1, start= start1,  data=data,
               control = nls.control(maxiter = 50 , warnOnly = TRUE),
               algorithm = "port")

  print("Halsey Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Halsey isotherm model

  errors <- function(y){
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))

# Graphical representation of the Halsey isotherm model

  ### Predicted parameter values
  parshal <- as.vector(coefficients(fit2))
  pars_nh <- parshal[1L];
  pars_Kh <- parshal[2L];

  rhs <- function(x){(exp((log(pars_Kh)-log(x))/pars_nh))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Halsey Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
