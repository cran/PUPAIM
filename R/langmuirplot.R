#' @title Langmuir Plot-deprecated
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the plot for the langmuir adsorption isotherm
#'
#' @name langmuirplot-deprecated
#' @usage langmuirplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{langmuirplot}:
#' For \code{langmuirplot}, use \code{\link{langmuiranalysis}}.
#'
#' @export
#' @export
langmuirplot <- function (Ce, Qe)
{ .Deprecated("langmuiranalysis")
  x <- 1/Ce
  y <- 1/Qe
  plot(x, y, main = "Langmuir Isotherm", xlab = "1 / Ce", ylab = "1 / Qe")
  fit86 <- lm(y ~ x)
  abline(fit86)
}
