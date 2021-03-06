% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harkinsjuraanalysis.R
\name{harkinsjuraanalysis}
\alias{harkinsjuraanalysis}
\title{HarkinsJura Isotherm Analysis Non-Linear Form}
\usage{
harkinsjuraanalysis(Ce, Qe)
}
\arguments{
\item{Ce}{the numerical value for the equilibrium capacity}

\item{Qe}{the numerical value for the adsorbed capacity}
}
\value{
the nonlinear regression and the parameters for the HarkinsJura isotherm
}
\description{
It assumes the possibility of multilayer adsorption on the surface of absorbents having heterogenous pore distribution (Ayawei, et al., 2017),(Gupta, et al., 2012)
}
\examples{
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
harkinsjuraanalysis(Ce, Qe)
}
\references{
Ayawei, N., Ebelegi, A. N., & Wankasi, D. (2017). Modelling and Interpretation of Adsorption Isotherms. Journal of Chemistry, 2017, 1-11. doi: 10.1155/2017/3039817

Gupta, V. K., Pathania, D., Agarwal, S., & Sharma, S. (2012). Removal of Cr(VI) onto Ficus carica biosorbent from water. Environmental Science and Pollution Research, 20(4), 2632-2644. doi:10.1007/s11356-012-1176-6
}
\author{
Raymond James L. Adame

C.C. Deocaris
}
