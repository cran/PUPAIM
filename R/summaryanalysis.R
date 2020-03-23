#' @title Summary of the Isotherm Analysis
#' @description 1-2 sentences Summarize the analysis for different isotherm models
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return summary of the nonlinear and linear fitting for different adsorption isotherm models
#' @importFrom stats "nls"
#' @importFrom utils "data"
#' @export

summaryanalysis <- function(Ce,Qe){
  fiveparamanalysis(Ce,Qe)
  bauduanalysis(Ce,Qe)
  BETanalysis(Ce,Qe)
  dubininradanalysis(Ce,Qe)
  fritzanalysis(Ce,Qe)
  harkinsjuraanalysis(Ce,Qe)
  henryanalysis(Ce,Qe)
  jossensanalysis(Ce,Qe)
  jovanovicanalysis(Ce,Qe)
  kahnanalysis(Ce,Qe)
  koblecarrigananalysis(Ce,Qe)
  langmuirFanalysis(Ce,Qe)
  mjanalysis(Ce,Qe)
  raudkepanalysis(Ce,Qe)
  redlichpanalysis(Ce,Qe)
  sipsanalysis(Ce,Qe)
  print(summary(summaryanalysis))
}
