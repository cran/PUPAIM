#' @title Marckzewski-Jaroniec Isotherm Nonlinear Analysis
#' @name marckzewskijaroniecanalysis
#' @description The Marczewski-Jaroniec Isotherm model has a resemblance to
#' Langmuir Isotherm model. It is developed on the basis of the supposition of
#' local Langmuir isotherm and adsorption energies distribution in the active
#' sites on adsorbent. This equation comprises all isotherm equations being an
#' extension of simple Langmuir Isotherm to single solute adsorption on heterogeneous
#' solids.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the absorbed capacity
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @return the nonlinear regression, parameters for Marckzewski-Jaroniec isotherm,
#' and model error analysis
#' @examples Qe <- c(0.19729, 0.34870, 0.61475, 0.74324, 0.88544, 0.89007, 0.91067, 0.91067, 0.96114)
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples marckzewskijaroniecanalysis(Ce,Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Marczewski, A. W., Derylo-Marczewska, A., and Jaroniec, M. (1986)
#' <doi:10.1016/0021-9797(86)90309-7M> Energetic heterogeneity and molecular
#' size effects in physical adsorption on solid surfaces. Journal of Colloid And
#' Interface Science, 109(2), 310-324.
#' @export

# Building the Marckzewski-Jaroniec isotherm nonlinear form
marckzewskijaroniecanalysis <- function (Ce, Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x,y)

# Marckzewski-Jaroniec isotherm nonlinear equation
  fit1 <- y ~ Qmax*((K * x)^n /(1+(K * x)^n))^(M/n)

# Setting of starting values

  N <- nrow(na.omit(data))
  K_max <- ((mean(y)/mean(x))+1)
  start <- data.frame(Qmax = seq(0, 1000, length.out = N),
                      K = seq(0, K_max, length.out = N),
                      M = seq(0, 1, length.out = N),
                      n = seq(0, 1, length.out = N))

# Fitting of the Marckzewski-Jaroniec isotherm via nls2

  suppressWarnings(fit2 <- nls2::nls2(fit1, start = start, data = data,
                 control = nls.control(maxiter = 1000,  warnOnly = T),
                 algorithm  = "default"))

  print("Marckzewski-Jaroniec Isotherm Nonlinear Analysis")
  print(summary(fit2))

  print("Aikake Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Marckzewski-Jaroniec isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Root Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))

# Graphical representation of the Marckzewski-Jaroniec isotherm model

  ### Predicted parameter values
  parsmarckj <- as.vector(coefficients(fit2))
  pars_Qmax <- parsmarckj[1L];
  pars_K <- parsmarckj[2L];
  pars_M <- parsmarckj[3L];
  pars_n <- parsmarckj[4L];

  rhs <- function(x){pars_Qmax*((pars_K*x)^pars_n /
                                  (1+(pars_K * x)^pars_n))^(pars_M/pars_n)}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Marckzewski-Jaroniec Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
