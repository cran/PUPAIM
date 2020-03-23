#' @title Elovich Isotherm Analysis Non-Linear Form
#' @description Elovich empirical adsorption model is based on the assumption of energetic heterogeneity of the adsorption sites. Moreover, elovich maximum adsorption capacity and Elovich constant can be calculated from the slopes and the intercepts of the plot ln(Qe/Ce) versus Qe.
#' @param Qe the numerical value for adsorbed capacity
#' @param Ce the numerical value for equilibrium concentration
#' @importFrom graphics "plot"
#' @importFrom graphics "abline"
#' @importFrom nls2 "nls2"

#' @importFrom Metrics "rmse" "mae" "mse" "rae"
#' @return the linear regression and the parameters for the Elovich isotherm
#' @examples Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
#' @examples Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
#' @examples elovichanalysis(Qe, Ce)
#' @author Edmundo B. De Guzman Jr.
#' @author C.C. Deocaris
#' @references Ayawei, N., Ebelegi, A.N., & Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 1-11. doi: 10.1155/2017/3039817
#' @references Gutierrez, L.G., Moreno-Pirajan, J.C., Guarin Romero, J.R. (2018). Kinetic and Equilibrium Study of the Adsorption of CO2 in Ultramicropores of Resorcinol-Formaldehyde Aerogels Obtained in Acidic and Basic Medium. Journal of Carbon Research, 2018, pp. 4. doi:10.3390/c4040052
#' @references Farouq, R., & Yousef, N.S. (2015). Equilibrium and Kinetics Studies of Adsorption of Copper(II) Ions on Natural Biosorbent. International Journal of Chemical Engineering and Applications, 2015, pp.332. DOI: 10.7763/IJCEA.2015.V6.503
#' @export
elovichanalysis <- function(Qe, Ce){
  x<-Qe
  y<-Qe/Ce
  data <- data.frame(x, y)
  mod207 <- (y~Ke*Qm*exp(x/Qm))
  start<- data.frame(Ke = c(0, 1000), Qm = c(0, 1000))
  set.seed(511)
  suppressWarnings(fit208<-nls2(mod207, start = start, control = nls.control(maxiter = 50, warnOnly = TRUE), algorithm = "port"))
  print("Elovich Non-Linear Isotherm Analysis")
  print(summary(fit208))
  N <- nrow(na.omit(data))
  error <- function(y){
    pv  <- (predict(fit208))
    rmse<- (rmse(y,predict(fit208)))
    mae <- (mae(y,predict(fit208)))
    mse <- (mse(y,predict(fit208)))
    rae <- (rae(y,predict(fit208)))
    PAIC <- AIC(fit208)
    PBIC <- BIC(fit208)
    SE <-(sqrt(sum(predict(fit208)-x)^2)/(N-2))
    colnames(y) <- rownames(y) <- colnames(y)
    list("Predicted Values"           = pv,
         "Relative Mean Square Error" = rmse,
         "Mean Absolute Error"        = mae,
         "Mean Squared Error"         = mse,
         "Relative Absolute Error"    = rae,
         "AIC"                        = PAIC,
         "BIC"                        = PBIC,
         "Standard Error Estimate"    = SE)
  }
  e <- error(y)
  print(e)
  plot(x,y, main = "Elovich Isotherm Model", xlab = "Qe",
       ylab = "Qe/Ce")
  lines(x, predict(fit208), lty=1,col="black",lwd=1)
  rsqq <- lm(y~predict(fit208))
  print(summary(rsqq))
}
