#' @title Brunauer-Emett-Teller (BET) Isotherm Linear Analysis
#' @name BET.LM
#' @description BET was particularly formulated to describe the multilayer
#' adsorption process in gas systems, but can also be employed to an aqueous
#' solution that relates the binding between layers because of the molecular
#' charge among them.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for BET isotherm, and model error analysis
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples BET.LM(Ce,Qe)
#' @author Jemimah Christine L. Mesias
#' @author Chester C. Deocaris
#' @references Brunauer, S., Emmett, P.H. and Teller, E. (1938) <doi:10.1021/ja01269a023>
#' Adsorption of Gases in Multimolecular Layers. Journal of the American Chemical Society, 60, 309-319.
#' @export

# Building the BET isotherm linear form
BET.LM <- function (Ce, Qe){

  x1 <- Ce
  y1 <- Qe
  data <- data.frame(x1, y1)

### BET isotherm nonlinear equation
  fit1 <- y1 ~ (CBET *x1)/((Cs-x1)*(1+((CBET-1)*(x1/Cs)))) ### Qmax is conditionally linear

### Setting of starting values
  start1 <- data.frame(CBET = seq(-50, 400, length.out = 50),
                       Cs = seq(-400, 50, length.out = 50))

### Fitting of the BET isotherm via nls2
  fit2 <- nls2::nls2(fit1, start = start1,  data=data,
               control = nls.control(maxiter = 45 , warnOnly = TRUE),
               algorithm = "plinear-random")

  param <- summary(fit2)
  BET_Cs <- param$coefficients[2]

# Establishing BET isotherm linear form
  x <- Ce
  y <- Ce/(Qe*(BET_Cs-Ce))

# Fitting of the BET isotherm linear form
  rhs <- function(Ce, CBET, Cs, Qmax) {
     1/Qmax + ((CBET-1)/(Qmax*CBET))*(Ce/Cs)
  }

  fit3 <- lm(y~x)

  print("BET Isotherm Linear Analysis")
  print(summary(fit3))

  ### y = a+bx
  c <- (summary(fit3))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

  ### Parameter values calculation
  Qmax <- 1/a
  print("Qmax")
  print(Qmax)

  CBET <-  1/(1-b*(Qmax*BET_Cs))
  print("CBET")
  print(CBET)

  Cs <- BET_Cs
  print("Cs")
  print(Cs)
# ---------------------------------
  print("Akaike Information Criterion")
  print(AIC(fit3))

  print("Bayesian Information Criterion")
  print(BIC(fit3))

# Error Analysis of the BET isotherm model

  errors <- function(y) {
    rmse <- Metrics::rmse(y, predict(fit3))
    mae <- Metrics::mae(y, predict(fit3))
    mse <- Metrics::mse(y, predict(fit3))
    rae <- Metrics::rae(y, predict(fit3))
    N <- nrow(na.omit(data))
    SE <- sqrt((sum(y-predict(fit3))^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Relative Mean Squared Error" = rmse,
         "Mean Absolute Error" = mae,
         "Mean Squared Error" = mse,
         "Relative Absolute Error" = rae,
         "Standard Error for the Regression S" = SE)
  }
  a <- errors(y)
  print(a)

  # Graphical representation of the Brunauer-Emett-Teller (BET) isotherm model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = "Ce",
         y = "Ce/Qe(Cs-Ce)",
         title = "Brunauer-Emett-Teller (BET) Isotherm Linear Model",
         caption = "PUPAIM") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}
