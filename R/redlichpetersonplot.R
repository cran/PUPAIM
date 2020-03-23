#' @title Redlich Peterson Isotherm Plot-deprecated
#' @description Plot of the analysis of Redlich Peterson Isotherm
#' @param Ce Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return Provides the plot for Redlich Peterson isotherm model
#'
#' @name rpplot-deprecated
#' @usage rpplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{rpplot}:
#' For \code{rpplot}, use \code{\link{redlichpanalysis}}.
#'
#' @export
rpplot<- function(Ce, Qe){
  .Deprecated("redlichpanalysis")
    x <- log10(Ce)
    y <- Qe
    fit75 <- lm(y~x)
    plot(y~x, xlab="Ce", ylab="Qe", main="Redlich-Peterson Analysis")
    abline(fit75, col="black")
}
