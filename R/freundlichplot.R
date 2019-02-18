#' Title
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the plot for the freundlich adsorption isotherm
#' @examples freundlichplot(c(1,2,3,4,5),c(1,2,3,4,5))
#' @export
#'
freundlichplot <- function (Ce, Qe)
{
  x <- log10(Ce)
  y <- log10(Ce)
  plot(x, y, main = "Feundlich Isotherm", xlab = "log(Ce)",
       ylab = "log(Qe)")
  fit85 <- lm(y ~ x)
  abline(fit85)
}

