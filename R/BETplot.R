#' @title BET Isotherm Plot-deprecated
#' @description Plot of the analysis of BET Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the plot for the BET adsorption isoth
#'
#' @name betplot-deprecated
#' @usage betplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{betplot}:
#' For \code{betplot}, use \code{\link{BETanalysis}}.
#'
#' @export
betplot <- function(Ce,Qe){
  x <- Ce
  y <- 1/Qe
  .Deprecated("BETanalysis")
  fit84 <- lm(y~x)
  plot(y~x, main="BET Isotherm Analysis", xlab= Ce, ylab = "1/Qe")
  abline(fit84)
}
