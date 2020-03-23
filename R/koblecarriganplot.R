#' @title Koble-Carrigan Isotherm Plot-deprecated
#' @description Plot of the analysis of Koble-Carrigan Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the plot for koble carrigan isotherm analysis
#'
#' @name kcarriganplot-deprecated
#' @usage kcarriganplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{kcarriganplot}:
#' For \code{kcarriganplot}, use \code{\link{koblecarrigananalysis}}.
#'
#' @export
kcarriganplot <- function(Ce, Qe){
    .Deprecated("koblecarrigananalysis")
    x <- 1/ Ce
    y <- 1/Qe
    plot(y~x, xlab="1/ce", ylab="1/qe", main="Koble-Carrigan analysis")
    fit61 <- lm(y~x)
     abline(fit61, col="black")
}
