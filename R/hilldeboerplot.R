#' @title Hill Deboer Isotherm-deprecated
#' @description Plot of the analysis of Hill Deboer Isotherm
#' @param theta a numeric vector consists of fractional coverage
#' @param Ce a numeric vector consists of equilibrium concentration
#' @return the LSRL plot for Hill Deboer isotherm
#'
#' @name hdplot-deprecated
#' @usage hdplot(theta, Ce)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{hdplot}:
#' For \code{hdplot}, use \code{\link{hilldeboeranalysis}}.
#'
#' @export
hdplot<- function(theta, Ce){
    .Deprecated("hilldeboeranalysis")
    x <- theta/(8.314*(273.15+25))
    y <- -log10(Ce)
    fit45 <- lm(y~x)
    plot(y~x, xlab="Theta", ylab="Ce", main="Hill-Deboer Analysis")
abline(fit45 , col="black")
}
