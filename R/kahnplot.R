#' @title Kahn Isotherm plot
#' @description Plot of the analysis of Kahn Isotherm
#' @param Ce the numerical value for the equilibrium capacity
#' @param Qe the numerical value for the adsorbed capacity
#'
#' @return the LSRL plot for Kahn isotherm analysis
#' @examples kahnplot(moringa$Ce,moringa$Qe)
#' @export
kahnplot <- function(Ce,Qe)
{
  x<-Ce
  y<-Ce/Qe
  plot(x, y, main = "Kahn Isotherm", xlab = "Ce", ylab= "Ce/Qe", type = "p")
  fit55 <- lm(y~x)
  abline(fit55, col = "green")
}

