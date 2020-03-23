#' @title Five Parameter Isotherm Plot-deprecated
#' @description Plot of the analysis of Five Parameter Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#' @importFrom graphics "plot"

#' @importFrom minpack.lm "nlsLM"
#' @return the linear regression and the parameters for the Five parameter isotherm analysis
#'
#' @name fivePplot-deprecated
#' @usage fivePplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{fivePplot}:
#' For \code{fivePplot}, use \code{\link{fiveparamanalysis}}.
#'
#' @export
fivePplot <- function(Ce,Qe){
  .Deprecated("fiveparamanalysis")
	x <- Ce^2
	y <- Qe
  fit7 <- lm(y~x)
	plot(y~x, main="Five parameter Analysis", xlab="Ce^2",ylab = "Qe")
	abline(fit7, col="black")
}
