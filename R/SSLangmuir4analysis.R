#' @title Langmuir Isotherm Nonlinear Analysis via selfStart and Langmuir Fourth Linear Model
#' @name SSLangmuir4analysis
#' @description The Langmuir isotherm is described to be the most useful and
#' simplest isotherm for both chemical adsorption and physical adsorption. It
#' assumes that there is uniform adsorption energy onto the monolayer surface
#' and that there would be no interaction between the adsorbate and the surface.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression via selfStart, initial starting values for parameters
#' based on Langmuir fourth linear model, predicted parameter values, and
#' model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples SSLangmuir4analysis(Ce,Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Langmuir, I. (1918) <doi:10.1021/ja01269a066> The adsorption of
#' gases on plane surfaces of glass, mics and platinum. Journal of the American
#' Chemical Society, 1361-1403.
#' @export

# Building the Langmuir isotherm nonlinear form
SSLangmuir4analysis <- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data1 <- data.frame(x, y)
  data2 <- data.frame(Ce,Qe, end = 1:nrow(data1))

  ### Initial starting values of the selfStart function
  parsinit4 <- getInitial(Qe ~ SSLangmuir4(Ce,Qmax,Kl), data = data2)

  ### Nonlinear fitting using selfStart function
  fit1 <- nls(Qe ~ SSLangmuir4(Ce,Qmax,Kl), data = data2 )

  ### Predicted paramaters from nonlinear fitting
  parsLangmuir4 <- as.vector(coefficients(fit1))
  pars_Qmax <- parsLangmuir4[1L];
  pars_Kl <- parsLangmuir4[2L];

  print("Initial Starting Value")
  print(parsinit4)

  print("Langmuir Isotherm Nonlinear Analysis")
  print(summary(fit1))

  print("Akaike Information Criterion")
  print(AIC(fit1))

  print("Bayesian Information Criterion")
  print(BIC(fit1))

  # Error analysis of the Langmuir isotherm model
  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae(y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data1))
    SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Root Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)
  
  rsqq <- lm(Qe~predict(fit1))
  print(summary(rsqq))

  # Graphical representation of the Langmuir isotherm model
  rhs <- function(x){((pars_Qmax*pars_Kl*x)/(1+(pars_Kl*x)))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data1, ggplot2::aes(x  = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Langmuir Isotherm Nonlinear Model",
         subtitle = "via selfStart based on fourth linear form",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5),
          plot.subtitle =ggplot2::element_text(hjust = 0.5))
}

