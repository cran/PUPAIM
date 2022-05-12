#' @title Radke-Prausnitz Isotherm Linear Analysis
#' @name radkepraustnitz.LM
#' @description The Radke-Prausnitz isotherm model has several important
#' properties which provides a good fit over a wide range of adsorbate
#' concentrations but more preferred in most adsorption systems at low adsorbate
#' concentration.
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @import nls2
#' @import Metrics
#' @import stats
#' @import ggplot2
#' @return the linear regression, parameters for Radke-Prausnitz isotherm,
#' and model error analysis
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples radkepraustnitz.LM(Ce,Qe)
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Radke, C. J. and Prausnitz, J. M. (1972) <doi:10.1021/i160044a003>
#' Adsorption of organic solutions from dilute aqueous solution on activated carbon,
#' Ind. Eng. Chem. Fund. 11 (1972) 445-451.
#' @export
#'

# Building the Radke-Praustnitz isotherm linear form
radkepraustnitz.LM <- function(Ce,Qe){

  x1 <- Ce
  y1 <- Qe
  data <- data.frame(x1, y1)

# Obtaining the model exponent
### Radke-Praustnitz isotherm nonlinear equation
  fit1 <- y1 ~ (Qmax*Krp*x1)/(1+(Krp*(x1^Nrp)))

### Setting of starting values
  start1 <- list(Qmax = 1, Krp = 1, Nrp = 1)

### Fitting of the Radke-Praustnitz isotherm via nls2

  fit2 <- nls2::nls2(fit1, start = start1, data=data,
                 control = nls.control(maxiter = 100 , warnOnly = TRUE),
                 algorithm = "port")

  param <- summary(fit2)
  expModel <- param$coefficients[3]

# Establishing Radke-Prausnitz isotherm linear form
  x <- Ce^expModel
  y <- Ce/Qe

  rhs <- function (x, Qmax, Krp) {
    1/(Qmax*Krp) + (1/Qmax)*x
    }

# Radke-Prausnitz isotherm linear fitting
  fit3 <- lm(y~x)

  print("Radke-Prausnitz Isotherm Linear Analysis")
  print(summary(fit3))

### y = a + bx
  c <- (summary(fit3))
  a <- c$coefficients[1]
  b <- c$coefficients[2]

### Parameter values calculation
  Qmax <- b^-1
  print("Qmax")
  print(Qmax)

  Krp <- a*Qmax
  print("Krp")
  print(Krp)

  Nrp <- expModel
  print("Nrp")
  print(Nrp)
# -------------------------------------
  AIC <- AIC(fit3)
  print("Aikake Information Criterion")
  print(AIC)

  BIC <- BIC(fit3)
  print("Bayesian Information Criterion")
  print(BIC)

# Error analysis of the Radke-Praustnitz isotherm linear model
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
  a<- errors(y)
  print(a)

# Graphical representation of the Sips isotherm linear model

  #### Plot details
  ggplot2::theme_set(ggplot2::theme_bw(10))
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(color ="#3498DB" ) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", se = F, color = "#D35400" ) +
    ggplot2::labs(x = expression(paste("Ce"^"Nrp")),
         y = "Ce/Qe",
         title = "Radke-Praustnitz Isotherm Linear Model",
         caption = "PUPAIM 0.3.0") +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5))
}

