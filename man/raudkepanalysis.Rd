% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raudkeprausniizanalysis.R
\name{raudkepanalysis}
\alias{raudkepanalysis}
\title{Radke-Prausniiz Isotherm Analysis Nonlinear Form}
\usage{
raudkepanalysis(Ce, Qe)
}
\arguments{
\item{Ce}{the numerical value for the equilibrium capacity}

\item{Qe}{the numerical value for the adsorbed capacity}
}
\value{
the nonlinear regression and the parameters for the Raudke-Prausniiz isotherm
}
\description{
Nonlinear form of the Radke-Prausniiz Equation, the original form. The Radke-Prausnitz isotherm model has several important properties which makes it more preferred in most adsorption systems at low adsorbate concentration.
}
\examples{
raudkepanalysis(moringa$Ce,moringa$Qe)
}
\references{
Newton, A., & Donbebe. (2017, September 5). Modelling and Interpretation of Adsorption Isotherms. Retrieved from https://doi.org/10.1155/2017/3039817

Khalid, A., Kazmi, M. et.al, (2015). Kinetic & Equilibrium Modelling of Copper Biosorption. Retrieved from http://journals.pu.edu.pk/journals/index.php/jfet/article/view/527
}
\author{
Princess Joyce DL Reyes

Neil Ross S. Alayon

C.C. Deocaris
}
