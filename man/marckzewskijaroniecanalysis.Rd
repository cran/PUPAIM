% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marckzewskijaroniecanalysis.r
\name{marckzewskijaroniecanalysis}
\alias{marckzewskijaroniecanalysis}
\title{Marckzewski-Jaroniec Isotherm Nonlinear Analysis}
\usage{
marckzewskijaroniecanalysis(Ce, Qe)
}
\arguments{
\item{Ce}{the numerical value for the equilibrium capacity}

\item{Qe}{the numerical value for the absorbed capacity}
}
\value{
the nonlinear regression, parameters for Marckzewski-Jaroniec isotherm,
and model error analysis
}
\description{
The Marczewski-Jaroniec Isotherm model has a resemblance to
Langmuir Isotherm model. It is developed on the basis of the supposition of
local Langmuir isotherm and adsorption energies distribution in the active
sites on adsorbent. This equation comprises all isotherm equations being an
extension of simple Langmuir Isotherm to single solute adsorption on heterogeneous
solids.
}
\examples{
Qe <- c(0.19729, 0.34870, 0.61475, 0.74324, 0.88544, 0.89007, 0.91067, 0.91067, 0.96114)
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
marckzewskijaroniecanalysis(Ce,Qe)
}
\references{
Marczewski, A. W., Derylo-Marczewska, A., and Jaroniec, M. (1986)
<doi:10.1016/0021-9797(86)90309-7M> Energetic heterogeneity and molecular
size effects in physical adsorption on solid surfaces. Journal of Colloid And
Interface Science, 109(2), 310-324.
}
\author{
Keith T. Ostan

Chester C. Deocaris
}
