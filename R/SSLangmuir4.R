#' @title selfStart using Langmuir Fourth Linear Model
#' @name SSLangmuir4
#' @description It calculates initial estimates for the model parameters from data
#' so nls has a greater chance of convergence.
#' @usage SSLangmuir4(Ce, Qmax,Kl)
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qmax the maximum adsorption capacity
#' @param Kl the numerical value for the adsorbed capacity
#' @import stats
#' @return initial starting values for parameters based on Langmuir fourth linear model
#' @author Keith T. Ostan
#' @author Chester C. Deocaris
#' @references Langmuir, I. (1918) <doi:10.1021/ja01269a066> The adsorption of
#' gases on plane surfaces of glass, mics and platinum. Journal of the American
#' Chemical Society, 1361-1403.
#' @references Chen, X. (2015) <doi:10.3390/info6010014> Modeling of Experimental
#' Adsorption Isotherm Data. 14–22.
#' @export

# Defining the selfStart function for Langmuir isotherm
SSLangmuir4 <- selfStart(

langmuir4<- function(Ce, Qmax, Kl) {
  (Qmax*Kl*Ce)/(1+(Kl*Ce))
  },

# Initial value routine using Langmuir linear form 4
  langmuir4.init <- function(mCall,LHS,data) {
    xy <- sortedXyData(mCall[["Ce"]], LHS, data)
    x <- xy[,"x"]
    y <- xy[,"y"]

    ### Linear regression on pseudo-x and pseudo-y
    pseudoY <- y/x
    pseudoX <- y
    coefs <- coef(lm(pseudoY ~ pseudoX))
    Kl <- -coefs[2]
    Qmax <- coefs[1]/coefs[2]
    value <- c(Qmax, Kl)
    names(value) <- mCall[c("Qmax","Kl")]
    value
  },
    parameters = c("Qmax","Kl"))



