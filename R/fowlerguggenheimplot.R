#' @title Fowler Guggenheim Isotherm Plot-deprecated
#' @description Plot of the analysis of Fowler Guggenheim Isotherm
#' @param theta the numerical value for the fractional coverage
#' @param Ce the numerical value for the equilibrium concentration
#' @importFrom graphics "abline"
#'@importFrom graphics "plot"
#'@importFrom stats "lm"
#'@importFrom minpack.lm "nlsLM"
#' @return the LSRL plot for Flory Huggins isotherm analysis
#'
#' @name fgplot-deprecated
#' @usage fgplot(theta, Ce)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{fgplot}:
#' For \code{fgplot}, use \code{\link{fowlerganalysis}}.
#'
#' @export
fgplot<- function(theta, Ce){
    .Deprecated("fowlerganalysis")
    x <- theta/24.46529
    y <- Ce
    fit25 <- lm(y~x)
    plot(y~x, xlab="Theta", ylab="Ce", main="Fowler-Guggenheim Analysis")
    abline(fit25,col="black")
}
