#' @title selfStart using Langmuir First Linear Model
#' @name SSLangmuir1
#' @description It calculates initial estimates for the model parameters from data
#' so nls has a greater chance of convergence.
#' @usage SSLangmuir1(Ce, Qmax,Kl)
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qmax the maximum adsorption capacity
#' @param Kl the numerical value for the adsorbed capacity
#' @import stats
#' @return initial starting values for parameters based on Langmuir first linear model
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Langmuir, I. (1918) <doi:10.1021/ja01269a066> The adsorption of
#' gases on plane surfaces of glass, mics and platinum. Journal of the American
#' Chemical Society, 1361-1403.
#' @references Chen, X. (2015) <doi:10.3390/info6010014> Modeling of Experimental
#' Adsorption Isotherm Data. 14–22.
#' @export

# Defining the selfStart function for Langmuir isotherm
SSLangmuir1 <- selfStart(

langmuir1<- function(Ce,Qmax,Kl) {
  (Qmax*Kl*Ce)/(1+(Kl*Ce))
  },

# Initial value routine using Langmuir linear form 1
  langmuir1.init <- function(mCall,LHS,data) {
    xy <- sortedXyData(mCall[["Ce"]], LHS, data)
    x <- xy[,"x"]
    y <- xy[,"y"]

    ### Linear regression on pseudo-x and pseudo-y
    pseudoY <- x/y
    pseudoX <- x
    coefs <- coef(lm(pseudoY ~ pseudoX))
    Qmax <- 1/coefs[2]
    Kl <- (1/Qmax)*(1/coefs[1])
    value <- c(Qmax, Kl)
    names(value) <- mCall[c("Qmax","Kl")]
    value
  },
    parameters = c("Qmax","Kl"))



