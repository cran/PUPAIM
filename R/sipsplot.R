#' @title Sips Isotherm Plot-deprecated
#' @description Plot of the analysis of Sips Isotherm
#' @param Ce the numerical value for the equilbrium concentration
#' @param Qe the numerical value for the adsorbed concentration
#' @return the LSRL of Sips analysis
#'
#' @name sipsplot-deprecated
#' @usage sipsplot(Qe, Ce)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL
#' @rdname PUPAIM-deprecated
#' @section \code{sipsplot}:
#' For \code{sipsplot}, use \code{\link{sipsanalysis}}.
#'
#' @export
sipsplot<- function(Qe,Ce){
  .Deprecated("sipsanalysis")
  x <- -log10(1/Qe)
	y <- log10(Ce)
	fit77 <- lm(y~x)
plot(y~x, xlab="ln (1/Qe)", ylab="ln Ce", main="Sips Isotherm Plot")
abline(fit77, col="black")
}
