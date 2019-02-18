#' @title Summary of Plots of the Isotherm Analysis
#' @description Summarize the plot for different isotherm models
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the summary of the LSRL plot for different adsorption isotherm
#' @export
summaryplots <- function(Ce,Qe){
  drplot(Ce,Qe)
  fivePplot(Ce,Qe)
  halseyplot(Ce,Qe)
  henryanalysisplot(Ce,Qe)
  hillplot(Ce,Qe)
  jossensplot(Ce,Qe)
  jovanovicplot(Ce,Qe)
  kahnplot(Ce,Qe)
  kcarriganplot(Ce,Qe)
  kiselevplot(Ce,Qe)
  marcJplot(Ce,Qe)
  radkePplot(Ce,Qe)
  rpplot(Ce,Qe)
  rpplot(Ce,Qe)
}
