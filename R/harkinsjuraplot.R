#' @title Harkins- Jura Isotherm Plot-deprecated
#' @description  Plot of the analysis of Harkins-Jura Isotherm
#' @param Ce the numerical value for a equilbrium concentration
#' @param Qe the numerical value for the adsorbed quantities
#' @return the LSRL plot for the Harkins Jura Isotherm
#'
#' @name hjplot-deprecated
#' @usage hjplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{hjplot}:
#' For \code{hjplot}, use \code{\link{harkinsjuraanalysis}}.
#'
#' @export
#' @export
hjplot <- function(Ce,Qe){
  .Deprecated("harkinsjuraanalysis")
	x4 <- log(Ce)
	y4 <- (1/(Qe)^2)
	plot(x4,y4, main="Harkin-Jura Isotherm Analysis" , xlab="log(ce)" , ylab="(1/(Qe)^2)" , type="p", col="blue", pch=10)
	fit37 <- lm(y4 ~ x4)
	abline(fit37, col="green")
}

