#' @title Weber Van Vliet Isotherm Plot-deprecated
#' @description Weber Van Vliet Isotherm Plot
#' @param Qe the numeric value for the adsorbed concentration
#' @param Ce the numeric value for the equilibrium concentration
#' @return the linear plot for the Weber Van Vliet Isotherm
#'
#' @name WVVplot-deprecated
#' @usage WVVplot(Qe, Ce)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{WVVplot}:
#' For \code{WVVplot}, use \code{\link{webervvanalysis}}.
#'
#' @export
WVVplot <- function(Qe, Ce){
  .Deprecated("webervvanalysis")
 x <- Qe
	y <- Ce
plot(y~x, main="WVV analysis", xlab="Ce", ylab="Qe")
fit83 <- lm(y~x)
abline(fit83, col="black")
}
