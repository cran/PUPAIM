#' @title Toth Isotherm Plot-deprecated
#' @description Plot of the analysis of Toth Isotherm
#' @param Ce the numeric value for the equilibrium concentration
#' @param theta the numeric value for the fractional coverage
#' @return the LSRL plot for a toth isotherm analysis
#'
#' @name tothplot-deprecated
#' @usage tothplot(Ce, theta)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{tothplot}:
#' For \code{tothplot}, use \code{\link{tothanalysis}}.
#'

#' @export
tothplot <- function(Ce,theta){
    .Deprecated("tothanalysis")
    x <- log10(Ce)
    y <- theta
    plot(x, y, xlab="ce", ylab="theta", main="Toth Isotherm Plot")
    fit80 <- lm(y~x)
     abline(fit80, col="black")
}
