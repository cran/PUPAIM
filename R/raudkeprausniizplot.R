#' @title Raudke Prausniiz Isotherm Plot-deprecated
#' @description Plot of the analysis of Raudke Prausniiz Isotherm
#' @param Ce Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return Provides the plot for Raudke Prausniiz isotherm model
#'
#' @name raudkePplot-deprecated
#' @usage raudkePplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{raudkePplot}:
#' For \code{raudkePplot}, use \code{\link{raudkepanalysis}}.
#'
#' @export
raudkePplot<- function(Ce,Qe){
  .Deprecated("raudkepanalysis")
	x <- 1/(1+Ce)^2
	y <- 1/Qe
	fit72<- lm(y~x)
plot(y~x, main="Raudke-Prausniiz Analysis", xlab="1/Ce", ylab="1/Qe")
abline(fit72,col="black")
}
