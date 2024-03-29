% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSLangmuir3.R
\name{SSLangmuir3}
\alias{SSLangmuir3}
\title{selfStart using Langmuir Third Linear Model}
\usage{
SSLangmuir3(Ce, Qmax,Kl)
}
\arguments{
\item{Ce}{the numerical value for the equilibrium capacity}

\item{Qmax}{the maximum adsorption capacity}

\item{Kl}{the numerical value for the adsorbed capacity}
}
\value{
initial starting values for parameters based on Langmuir third linear model
}
\description{
It calculates initial estimates for the model parameters from data
so nls has a greater chance of convergence.
}
\references{
Langmuir, I. (1918) <doi:10.1021/ja01269a066> The adsorption of
gases on plane surfaces of glass, mics and platinum. Journal of the American
Chemical Society, 1361-1403.

Chen, X. (2015) <doi:10.3390/info6010014> Modeling of Experimental
Adsorption Isotherm Data. 14–22.
}
\author{
Jemimah Christine L. Mesias

Chester C. Deocaris
}
