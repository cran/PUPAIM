#'@title Redlich-Peterson Isotherm Nonlinear Analysis
#'@name redlichpetersonanalysis
#'@description Redlich-Peterson isotherm model has an exponential function which
#'can be found in the denominator and in the numerator, it has a linear dependence
#'on the concentration denoting the adsorption equilibrium depending on a wide
#'range of concentration
#'@param Ce the numerical value for the equilibrium capacity
#'@param Qe the numerical value for the adsorbed capacity
#'@import nls2
#'@import Metrics
#'@import stats
#'@import ggplot2
#'@return the nonlinear regression, parameters for Redlich-Peterson isotherm, and model error analysis
#'@examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#'@examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#'@examples redlichpetersonanalysis(Ce,Qe)
#'@author Keith T. Ostan
#'@author Chester C. Deocaris
#'@references Peterson, D. L. and Redlich, O.(1959) <doi:10.1021/j150576a611> A useful adsorption isotherm.
#'J PhysChem US;63(6):1024. Research, vol. 6, no. 1, pp. 265-276, 2012.
#'@export
#'

# Building the Redlich-Peterson isotherm nonlinear form
redlichpetersonanalysis <- function(Ce, Qe){

  x <- Ce
  y <- Qe
  data <- data.frame(x, y)

# Redlich-Peterson isotherm nonlinear equation
  fit1 <- y ~ (Krp*x)/(1+(Arp*(x^g)))

# Setting of starting values
  start1 <- list(Krp = 1, Arp = 1, g = 1)

# Fitting of the Redlich-Peterson isotherm via nls2

  fit2 <- nls2::nls2(fit1,start = start1, data=data,
                 control = nls.control(maxiter = 50, warnOnly = TRUE),
                 algorithm = "port")

  print("Redlich-Peterson Isotherm Nonlinear Analysis")
  print(summary(fit2))

  AIC <- AIC(fit2)
  print("Aikake Information Criterion")
  print(AIC)

  BIC <- BIC(fit2)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Redlich-Peterson isotherm model

  errors <- function(y){
    rmse <-Metrics::rmse(y, predict(fit2))
    mae <- Metrics::mae(y, predict(fit2))
    mse <- Metrics::mse(y, predict(fit2))
    rae <- Metrics::rae(y, predict(fit2))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit2))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Squared Error" = rmse,
           "Mean Absolute Error" = mae,
           "Mean Squared Error" = mse,
           "Relative Absolute error" = rae,
           "Standard Error for the Regression S" = SE)
    }
    a <- errors(y)
    print(a)
    
    rsqq <- lm(Qe~predict(fit2))
    print(summary(rsqq))

# Graphical representation of the Redlich-Peterson isotherm model

    ### Predicted parameter values
    parsredlichP <- as.vector(coefficients(fit2))
    pars_Krp <- parsredlichP[1L];
    pars_Arp <- parsredlichP[2L];
    pars_g <- parsredlichP[3L]

    rhs <- function(x){((pars_Krp*x)/(1+(pars_Arp*(x^pars_g))))}

    #### Plot details
    ggplot2::theme_set(ggplot2::theme_bw(10))
    ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
      ggplot2::geom_function(color = "#D35400", fun = rhs ) +
      ggplot2::labs(x = "Ce",
           y = "Qe",
           title = "Redlich-Peterson Isotherm Nonlinear Model",
           caption = "PUPAIM") +
      ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
