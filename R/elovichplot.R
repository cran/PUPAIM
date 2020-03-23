#' @title Elovich Isotherm Plot-deprecated
#' @description Plot of the analysis of Elovich Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#' @importFrom graphics "plot"
#' @importFrom minpack.lm "nlsLM"
#' @return the LSRL plot for Elovich isotherm analysis
#'
#' @name elovichplot-deprecated
#' @usage elovichplot(Qe, Ce)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{elovichplot}:
#' For \code{elovichplot}, use \code{\link{elovichanalysis}}.
#'
#' @export
elovichplot <- function(Qe,Ce){
  .Deprecated("elovichanalysis")
  x <- Qe
	y <- log10(Qe/Ce)
	plot(x, y, main ="Elovich Isotherm Plot", xlab= "Ce", ylab= "log10(Qe/Ce)", type="p", col="blue", pch=10)
	fit19 <- lm(y~x)
	abline(fit19, col="black")
}
