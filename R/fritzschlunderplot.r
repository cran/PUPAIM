#' @title Fritz Schlunder Isotherm Plot-deprecated
#' @description Plot of the analysis of Fritz Schlunder Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the LSRL plot fot Fritz Schlunder isotherm analysis
#'
#' @name fsplot-deprecated
#' @usage fsplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{fsplot}:
#' For \code{fsplot}, use \code{\link{fritzanalysis}}.
#'
#' @export
fsplot <- function(Ce,Qe){
  .Deprecated("fritzanalysis")
  x <- Ce
	y <- Qe
	plot(x, y, main ="Fritz-Schlunder Isotherm Plot", xlab= "Ce", ylab= "Qe", type="p", col="blue", pch=5)
	fit31 <- lm(y~x)
	abline(fit31, col="black")
}
