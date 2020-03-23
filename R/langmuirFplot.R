#' @title Langmuir-Freundlich Isotherm Plot-deprecated
#' @description Plot of the analysis of Langmuir-FreundLich Isotherm
#' @param Ce the numerical value for the equilbrium concentration
#' @param Qe the numerical value for the adsorbed concentration
#' @return the LSRL for the Langmuir-Freundlich isotherm analysis
#'
#' @name langmuirFplot-deprecated
#' @usage langmuirFplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{langmuirFplot}:
#' For \code{langmuirFplot}, use \code{\link{langmuirFanalysis}}.
#'
#' @export
langmuirFplot <- function(Ce,Qe){
  .Deprecated("langmuirFanalysis")
  x <- 1/(Ce)^2
	y <- Qe
	fit66 <- lm(y~x)
	plot(y~x, main="Langmuir-Freundlich Analysis", xlab="1/Ce^2", ylab = "Qe")
	abline(fit66,col="black")
}
