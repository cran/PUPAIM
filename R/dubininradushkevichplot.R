#' @title Dubinin-Radushkevich Isotherm Plot-deprecated
#' @description Plot of the analysis of Dubinun-Radushkevich Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @importFrom graphics "abline"
#' @importFrom graphics "plot"
#' @importFrom stats "lm"
#' @importFrom minpack.lm "nlsLM"
#' @return the LSRL plot of Dubinin-Radushkevich isotherm analysis
#'
#' @name drplot-deprecated
#' @usage drplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{drplot}:
#' For \code{drplot}, use \code{\link{dubininradanalysis}}.
#'
#' @export
drplot<- function(Ce, Qe){
    .Deprecated("dubininradanalysis")
    x <- Ce^2
    y <- log10(Qe)
    fit15<- lm(y~x)
    plot(y~x, xlab="Ce^2", ylab="Qe", main="Dubinin-Radushkevich Analysis")
   abline(fit15, col="black")
}
