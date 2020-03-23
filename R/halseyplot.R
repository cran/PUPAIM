#' @title Halsey Isotherm Plot-deprecated
#' @description Plot of the analysis of Halsey Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the LSRL plot for Halsey isotherm analysis
#'
#' @name halseyplot-deprecated
#' @usage halseyplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{halseyplot}:
#' For \code{halseyplot}, use \code{\link{halseyanalysis}}.
#'
#' @export
halseyplot <- function(Ce,Qe){
    .Deprecated("halseyanalysis")
    x<-log10(Ce)
    y <- Qe
    plot(x, y, xlab="ln ce", ylab="Qe", main="Halsey Analysis")
    fit34<- lm(y~x)
     abline(fit34, col="black")
}
