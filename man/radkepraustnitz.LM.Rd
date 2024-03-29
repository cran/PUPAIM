% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radkepraustnitzLM.R
\name{radkepraustnitz.LM}
\alias{radkepraustnitz.LM}
\title{Radke-Prausnitz Isotherm Linear Analysis}
\usage{
radkepraustnitz.LM(Ce, Qe)
}
\arguments{
\item{Ce}{the numerical value for the equilibrium capacity}

\item{Qe}{the numerical value for the adsorbed capacity}
}
\value{
the linear regression, parameters for Radke-Prausnitz isotherm,
and model error analysis
}
\description{
The Radke-Prausnitz isotherm model has several important
properties which provides a good fit over a wide range of adsorbate
concentrations but more preferred in most adsorption systems at low adsorbate
concentration.
}
\examples{
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
radkepraustnitz.LM(Ce,Qe)
}
\references{
Radke, C. J. and Prausnitz, J. M. (1972) <doi:10.1021/i160044a003>
Adsorption of organic solutions from dilute aqueous solution on activated carbon,
Ind. Eng. Chem. Fund. 11 (1972) 445-451.
}
\author{
Keith T. Ostan

Chester C. Deocaris
}
