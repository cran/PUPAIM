## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PUPAIM)

## ----Adsorption Plot, warning=FALSE-------------------------------------------
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
data <- data.frame(Ce,Qe)

plot(data,main = "Adsorption Data", xlab = "Ce", 
         ylab = "Qe")

## ----One-Parameter Adsorption Isotherm Model, warning=FALSE-------------------
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
data <- data.frame(Ce,Qe)

henry.LM(Ce, Qe)

## ----Two-Parameter Adsorption Isotherm Model, warning=FALSE-------------------
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
data <- data.frame(Ce,Qe)

hill.LM(Ce, Qe)
hillanalysis(Ce, Qe)

## ----Three-Parameter Adsorption Isotherm Model, warning=FALSE-----------------
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
data <- data.frame(Ce,Qe)

redlichpeterson.LM(Qe, Ce)
redlichpetersonanalysis(Qe, Ce)

## ----Four-Parameter Adsorption Isotherm Model, warning=FALSE------------------
Ce <- c(0.01353, 0.04648, 0.13239, 0.27714, 0.41600, 0.63607, 0.80435, 1.10327, 1.58223)
Qe <- c(0.03409, 0.06025, 0.10622, 0.12842, 0.15299, 0.15379, 0.15735, 0.15735, 0.16607)
data <- data.frame(Ce,Qe)

FS4analysis(Ce, Qe)

