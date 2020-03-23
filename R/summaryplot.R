#' @title Summary of Plots of the Isotherm Analysis-deprecated
#' @description Summarize the plot for different isotherm models
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the summary of the LSRL plot for different adsorption isotherm
#'
#' @name summaryplots-deprecated
#' @usage summaryplots(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL
#' @rdname PUPAIM-deprecated
#' @section \code{summaryplots}:
#' For \code{summaryplots}, use \code{\link{summaryanalysis}}.
#'
#' @export
summaryplots <- function(Ce,Qe){
  .Deprecated("summaryanalysis")
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
  rpplot(Ce,Qe)
}
