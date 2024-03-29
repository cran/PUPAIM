% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kahnanalysis.R
\name{kahnanalysis}
\alias{kahnanalysis}
\title{Kahn Isotherm Non-Linear Analysis}
\usage{
kahnanalysis(Ce, Qe)
}
\arguments{
\item{Ce}{the numerical value for the equilibrium capacity}

\item{Qe}{the numerical value the absorbed capacity}
}
\value{
the nonlinear regression, parameters for the Kahn isotherm, and
model error analysis
}
\description{
A generalized model recommended for pure solutions, in which
both extremes, Langmuir and Freundlich, can be represented. This isotherm was
developed to cater to both the single- and multi-component adsorption systems.
}
\examples{
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
kahnanalysis(Ce, Qe)
}
\references{
Khan, A. R., Al-Waheab, I. R., and Al-Haddad, A. (1996) <doi:10.1080/09593331708616356> A
generalized equation for adsorption isotherms for multi-component organic
pollutants in dilute aqueous solution. Environmental Technology (United Kingdom), 17(1), 13-23.
}
\author{
Paul Angelo C. Manlapaz

Chester C. Deocaris
}
