#' @title Jossens Isotherm Plot-deprecated
#' @description Plot of the analysis of Jossens Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#' @return the LSRL plot for Jossen's isotherm
#'
#' @name jossensplot-deprecated
#' @usage jossensplot(Ce, Qe)
#' @seealso \code{\link{PUPAIM-deprecated}}
#' @keywords internal
NULL

#' @rdname PUPAIM-deprecated
#' @section \code{jossensplot}:
#' For \code{jossensplot}, use \code{\link{jossensanalysis}}.
#'
#' @export
jossensplot<- function(Ce, Qe){
.Deprecated("jossensanalysis")
      x <- (Qe)
    y <- log10(Ce/Qe)
    fit49 <- lm(y~x)
    plot(y~x, xlab="Qe", ylab="ln(Ce/Qe)", main="Jossens Analysis")
  abline(fit49, col="blue")
   	}


