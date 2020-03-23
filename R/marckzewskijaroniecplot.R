#' @title Marckzewski Jaroniec Isotherm Plot-deprecated
#' @description Plot of the analysis of Marckzewski Jaroniec Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the plot for the Marckzewski Jaroniec isotherm analysis
#'
#' @name marcJplot-deprecated
#' @usage marcJplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{marcJplot}:
#' For \code{marcJplot}, use \code{\link{mjanalysis}}.
#'
#' @export
marcJplot <- function(Ce, Qe){
  .Deprecated("mjanalysis")
x <- (1/Ce^2)
y <- Qe
plot(y~x)
fit69 <- lm(y~x)
abline(fit69, main="M-J plot", xlab="1/Ce", ylab="Qe")
}
