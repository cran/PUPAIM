#' @title selfStart using Langmuir Second Linear Model
#' @name SSLangmuir2
#' @description It calculates initial estimates for the model parameters from data
#' so nls has a greater chance of convergence.
#' @usage SSLangmuir2(Ce, Qmax,Kl)
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qmax the maximum adsorption capacity
#' @param Kl the numerical value for the adsorbed capacity
#' @import stats
#' @return initial starting values for parameters based on Langmuir second linear model
#' @author Paul Angelo C. Manlapaz
#' @author Chester C. Deocaris
#' @references Langmuir, I. (1918) <doi:10.1021/ja01269a066> The adsorption of
#' gases on plane surfaces of glass, mics and platinum. Journal of the American
#' Chemical Society, 1361-1403.
#' @references Chen, X. (2015) <doi:10.3390/info6010014> Modeling of Experimental
#' Adsorption Isotherm Data. 14–22.
#' @export


# Defining the selfStart function for Langmuir isotherm
SSLangmuir2 <- selfStart(

langmuir2<- function(Ce, Qmax , Kl) {
  (Qmax*Kl*Ce)/(1+(Kl*Ce))
  },

# Initial value routine using Langmuir linear form 2
  langmuir2.init <- function(mCall,LHS,data) {
    xy <- sortedXyData(mCall[["Ce"]], LHS, data)
    x <- xy[,"x"]
    y <- xy[,"y"]

    ### Linear regression on pseudo-x and pseudo-y
    pseudoY <- 1/y
    pseudoX <- 1/x
    coefs <- coef(lm(pseudoY ~ pseudoX))
    Qmax <- 1/coefs[1]
    Kl <- (1/coefs[2])*1/Qmax
    value <- c(Qmax, Kl)
    names(value) <- mCall[c("Qmax","Kl")]
    value
  },
    parameters = c("Qmax","Kl"))



