#' @title Jovanovic Isotherm Plot-deprecated
#' @description Plot of the analysis of Jovanovic Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity

#' @return the LSRL plot for Jovanovic isotherm analysis
#'
#' @name jovanovicplot-deprecated
#' @usage jovanovicplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{jovanovicplot}:
#' For \code{jovanovicplot}, use \code{\link{jovanovicanalysis}}.
#'
#' @export
jovanovicplot <- function(Ce,Qe)
  {
  .Deprecated("jovanovicanalysis")
	x5<-Ce
	y5<-log10(Qe)
	plot(x5, y5, main = "Jovanovic Isotherm", xlab = "Ce", ylab= "ln(Qe)", type = "p", col ="blue", pch =10)
	fit52 <- lm(y5~x5)
	abline(fit52, col = "green")
  }

