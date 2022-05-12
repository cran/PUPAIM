#' @title Flory-Huggins Isotherm Non-Linear Analysis
#' @name floryhugginsanalysis
#' @description Flory-Huggins isotherm model describes the degree of surface coverage
#' characteristics of the adsorbate on the adsorbent. It describes the nature of the
#' adsorption process regarding the feasibility and spontaneity of the process. The theory
#' of the Flory-Huggins provides the mathematical model for the polymer blends'
#' thermodynamics.
#' @param Ce is equal to Co which is the numeric value for the initial concentration
#' @param theta is the fractional surface coverage
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Flory-Huggins isotherm, and model
#' error analysis
#' @examples theta <- c(0.19729, 0.34870, 0.61475, 0.74324, 0.88544, 0.89007, 0.91067, 0.91067, 0.96114)
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples floryhugginsanalysis(theta, Ce)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Flory, P. J. (1971). Principles of polymer chemistry. Cornell Univ.Pr.
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export
#'

# Building the Flory-Huggins isotherm nonlinear form
floryhugginsanalysis <- function(theta, Ce){

  x <- theta
  y <- 1/Ce
  data <- data.frame(x, y)

# Flory-Huggins isotherm nonlinear equation
  fit1 <- y ~ (1/x)*(KFH*(1-x)^nFH)

# Setting of starting values
  start1 <- list(KFH = 100 , nFH = 1)

# Fitting of the Flory-Huggins isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1,  data=data,
                 control= nls.control(maxiter= 100, warnOnly=TRUE),
                 algorithm= "port")

  print("Flory HUggins Parameters")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Flory-Huggins isotherm model

  errors <- function(y) {
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

  # Graphical representation of the Flory-Huggins isotherm model

  ### Predicted parameter values
  parsFloryHuggins <- as.vector(coefficients(fit2))
  pars_KFH <- parsFloryHuggins[1L];
  pars_nFH <- parsFloryHuggins[2L];

  rhs <- function(x) ((1/x)*(pars_KFH *(1-x)^pars_nFH))

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = expression(paste(theta)),
         y = "1/Ce",
         title = "Flory-Huggins Isotherm Nonlinear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
