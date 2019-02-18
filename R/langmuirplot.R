#' Langmuir Plot
#'
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the plot for the langmuir adsorption isotherm
#' @examples langmuirplot(moringa$Ce,moringa$Qe)
#' @export
langmuirplot <- function (Ce, Qe)
{
  x <- 1/Ce
  y <- 1/Qe
  plot(x, y, main = "Langmuir Isotherm", xlab = "1 / Ce", ylab = "1 / Qe")
  fit86 <- lm(y ~ x)
  abline(fit86)
}
