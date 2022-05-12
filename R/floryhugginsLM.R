#' @title Flory-Huggins Isotherm Linear Analysis
#' @name floryhuggins.LM
#' @description Flory-Huggins isotherm model describes the degree of surface
#' coverage characteristics of the adsorbate on the adsorbent. It describes the
#' nature of the adsorption process regarding the feasibility and spontaneity of
#' the process. The theory of the Flory-Huggins provides the mathematical model
#' for the polymer blends' thermodynamics.
#' @param Ce the numerical value for the equilibrium capacity
#' @param theta is theta fractional surface coverage
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Flory-Huggins isotherm, and
#' model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607,
#' 0.80435, 1.10327, 1.58223)
#' @examples theta  <- c(0.1972984, 0.3487013, 0.6147560, 0.7432401, 0.8854408,
#' 0.8900708, 0.9106746, 0.9106746, 0.9611422)
#' @examples floryhuggins.LM (theta,Ce)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Flory, P. J. (1971). Principles of polymer chemistry. Cornell Univ.Pr.
#' @references Foo, K. Y., &amp; Hameed, B. H. (2009, September 13).
#' <doi:10.1016/j.cej.2009.09.013> Insights into the modeling of adsorption isotherm
#' systems. Chemical Engineering Journal.
#' @export

# Building the Flory-Huggins isotherm linear form
floryhuggins.LM <- function(theta,Ce){

  x <- 1 - theta
  y <- log10(theta/Ce)
  data <- data.frame(x, y)

# Fitting of the Flory-Huggins isotherm linear form
rhs <- function(x, KFH, nFH) {
    log(KFH)+(nFH)*log(1-theta)
  }

  fit1 <- lm(y~x)

  print("Flory-Huggins Isotherm Analysis")
  print(summary(fit1))

  ### y = a + bx
  c <- summary(fit1)
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  KFH <- 10^(a)
  print("KFH")
  print(KFH)

  nFH <- b
  print("nFH")
  print(nFH)

# ---------------------------------
  AIC <- AIC(fit1)
  print("Akaike Information Criterion")
  print(AIC)

  BIC <- BIC(fit1)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Flory-Huggins isotherm model
  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit1))
    mae <- Metrics::mae(y, predict(fit1))
    mse <- Metrics::mse(y, predict(fit1))
    rae <- Metrics::rae(y, predict(fit1))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit1))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

# Graphical representation of the Flory-Huggins isotherm linear model

#### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste("1-", theta)),
         y = expression(paste("log(", theta,"/Ce)")),
         title = "Flory-Huggins Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
