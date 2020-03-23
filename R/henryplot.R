#' @title Henry Isotherm Plot-deprecated
#' @description Plot of the analysis of Henry Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the LSRL plot for Henry isotherm analysis
#'
#' @name henryanalysisplot-deprecated
#' @usage henryanalysisplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{henryanalysisplot}:
#' For \code{henryanalysisplot}, use \code{\link{henryanalysis}}.
#'
#' @export
henryanalysisplot <- function(Ce,Qe){
  .Deprecated("henryanalysis")
	x1 <- Ce
	y1 <- Qe
	plot(x1, y1, main="Henry Isotherm Analysis" , xlab= "Ce" , ylab="Qe", type="p" , col="blue", pch=10)
	fit40 <- lm(y1~x1)
	abline(fit40, col="green")
}
