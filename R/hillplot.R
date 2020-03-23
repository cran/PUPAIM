#' @title Hill Isotherm Plot-deprecated
#' @description Plot of the analysis of Hill Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the LSRL plot for Hill isotherm analysis
#'
#' @name hillplot-deprecated
#' @usage hillplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{hillplot}:
#' For \code{hillplot}, use \code{\link{hillanalysis}}.
#'
#' @export
hillplot<- function(Ce,Qe){
    .Deprecated("hillanalysis")
    x <- log(Ce)
    y <- log(Qe/.5-Qe)
    fit46 <- lm(y~x)
    plot(y~x, xlab="log(Ce)", ylab="Qe", main="Hill Analysis")
    abline(fit46,col="black")
}
