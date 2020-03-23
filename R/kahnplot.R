#' @title Kahn Isotherm plot-deprecated
#' @description Plot of the analysis of Kahn Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the LSRL plot for Kahn isotherm analysis
#'
#' @name kahnplot-deprecated
#' @usage kahnplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{kahnplot}:
#' For \code{kahnplot}, use \code{\link{kahnanalysis}}.
#'
#' @export
kahnplot <- function(Ce,Qe)
{
  .Deprecated("kahnanalysis")
  x<-Ce
  y<-Ce/Qe
  plot(x, y, main = "Kahn Isotherm", xlab = "Ce", ylab= "Ce/Qe", type = "p")
  fit55 <- lm(y~x)
  abline(fit55, col = "green")
}

