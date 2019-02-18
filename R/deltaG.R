#' @title Gibbs free energy of Adsorption
#' @description Defines the spontaneity of an adsoprtion process base from the value. If it is negative is is spontaneous and if the result is positive,the reaction is not spontaneous.
#'@param t temperature used in the experiment
#'@param K equilirium constant for the adsorption process
#'
#' @return the Gibbs free energy value for the given adsorption process in terms of Kj(kilojoule) per mole
#' @examples deltaG(25,1234)
#' @export
deltaG <- function(t,K){
 x <- 8.314*(273.15+t)
  y <- log10(K)
(-(x*y))/1000
}
