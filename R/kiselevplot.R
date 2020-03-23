#' @title Kiselev Isotherm Plot-deprecated
#' @description Plot of the analysis of Kiselev Isotherm
#' @param Kh the numerical value for the kiselev constant
#' @param Ce the numerical value for the equilibrium capacity
#' @return the plot for kiselev isotherm analysis
#'
#' @name kiselevplot-deprecated
#' @usage kiselevplot(Qe, Ce)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{kiselevplot}:
#' For \code{kiselevplot}, use \code{\link{kiselevanalysis}}.
#'
#' @export
kiselevplot <- function(Qe,Ce){
  .Deprecated("kiselevanalysis")
  x<- Qe
  y<-1/Ce
     plot(x, y, xlab="Qe", ylab="1/Ce", main="Kiselev Analysis")
    fit58 <- lm(y~x)
     abline(fit58, col="black")
}
