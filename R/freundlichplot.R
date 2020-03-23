#' @title Freundlich Isotherm plot-deprecated
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the plot for the freundlich adsorption isotherm
#'
#' @name freundlichplot-deprecated
#' @usage freundlichplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{freundlichplot}:
#' For \code{freundlichplot}, use \code{\link{freundlichanalysis}}.
#'
#' @export
freundlichplot <- function (Ce, Qe)
{ .Deprecated("freundlichanalysis")
  x <- log10(Ce)
  y <- log10(Ce)
  plot(x, y, main = "Freundlich Isotherm", xlab = "log(Ce)",
       ylab = "log(Qe)")
  fit85 <- lm(y ~ x)
  abline(fit85)
}

