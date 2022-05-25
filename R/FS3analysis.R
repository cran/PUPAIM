#' @title Fritz-Schlunder Three Parameter Non-Linear Analysis
#' @name FS3analysis
#' @description The Fritz-Schlunder isotherm model is an empirical expression that
#' can fit over an extensive range of experimental results as a result of the huge
#' number of coefficients in their adsorption isotherm.
#' @param Ce the numerical value for equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the nonlinear regression, parameters for Fritz-Schlunder three Parameter
#' isotherm, and model error analysis
#' @examples Ce <- c(0.9613, 1.0895, 1.5378, 1.9862, 3.3314, 7.8153, 11.4024, 15.8862)
#' @examples Qe <- c(2.5546, 4.4150, 5.8558, 7.1387, 8.8092, 13.1921, 15.7966, 18.4483)
#' @examples FS3analysis(Ce,Qe)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Fritz, W., and Schluender, E. U. (1974) <doi:z10.1016/0009-2509(74)80128-4> Simultaneous adsorption
#' equilibria of organic solutes in dilute aqueous solutions on activated carbon.
#' Chemical Engineering Science, 29(5), 1279-1282.
#' @export

# Building the Fritz-Schlunder Three Parameter isotherm nonlinear form
FS3analysis <- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# Fritz-Schlunder Three Parameter isotherm nonlinear equation
  fit1 <- y ~ ((QmaxFS * KFS * x)/(1 + QmaxFS * x ^ (MFS)))

# Setting of starting values
  start1 <- data.frame(QmaxFS = 1, KFS = 1, MFS = 1)

# Fitting of Fritz-Schlunder Three Parameter isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter=10000 , warnOnly = TRUE),
                 algorithm = "port")

  print("Fritz-Schlunder Three Parameter Isotherm Non-linear Analysis")
  print(summary(fit2))

  print("Akaike Information Criterion")
  print(AIC(fit2))

  print("Bayesian Information Criterion")
  print(BIC(fit2))

# Error analysis of the Fritz-Schlunder Three Parameter isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Root Mean Squared Error" = rmse, "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Root Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  rsqq <- lm(Qe~predict(fit2))
  print(summary(rsqq))

  # Graphical representation of the Fritz-Schlunder Three Parameter isotherm model

  ### Predicted parameter values
  parsfritzthree <- as.vector(coefficients(fit2))
  pars_QmaxFS <- parsfritzthree[1L];
  pars_KFS <- parsfritzthree[2L];
  pars_MFS <- parsfritzthree[3L];

  rhs <- function (x){((pars_QmaxFS*pars_KFS*x)/(1+pars_QmaxFS*x^(pars_MFS)))}

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_function(color = "#D35400", fun = rhs ) +
    ggplot2::labs(x = "Ce",
         y = "Qe",
         title = "Fritz-Schlunder (III) Isotherm Nonlinear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

