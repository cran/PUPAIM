% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BETLM.r
\name{BET.LM}
\alias{BET.LM}
\title{Brunauer-Emett-Teller (BET) Isotherm Linear Analysis}
\usage{
BET.LM(Ce, Qe)
}
\arguments{
\item{Ce}{the numerical value for the equilibrium capacity}

\item{Qe}{the numerical value for the adsorbed capacity}
}
\value{
the linear regression, parameters for BET isotherm, and model error analysis
}
\description{
BET was particularly formulated to describe the multilayer
adsorption process in gas systems, but can also be employed to an aqueous
solution that relates the binding between layers because of the molecular
charge among them.
}
\examples{
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
BET.LM(Ce,Qe)
}
\references{
Brunauer, S., Emmett, P.H. and Teller, E. (1938) <doi:10.1021/ja01269a023>
Adsorption of Gases in Multimolecular Layers. Journal of the American Chemical Society, 60, 309-319.
}
\author{
Jemimah Christine L. Mesias

Chester C. Deocaris
}
