% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/halseyanalysis.R
\name{halseyanalysis}
\alias{halseyanalysis}
\title{Halsey Isotherm Non-Linear Analysis}
\usage{
halseyanalysis(Ce, Qe)
}
\arguments{
\item{Ce}{the numerical value for the equilibrium capacity}

\item{Qe}{the numerical value for the adsorbed capacity}
}
\value{
the nonlinear regression, parameters for the Halsey isotherm, and
model error analysis
}
\description{
A multilayer adsorption isotherm model which is suited for
adsorption of adsorbate ions at a distance that is relatively large from the
surface.
}
\examples{
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
halseyanalysis(Ce, Qe)
}
\references{
Halsey, G., and Taylor, H. S. (1947) <doi:10.1063/1.1746618> The adsorption of
hydrogen on tungsten powders. The Journal of Chemical Physics, 15(9), 624-630.
}
\author{
Paul Angelo C. Manlapaz

Chester C. Deocaris
}
