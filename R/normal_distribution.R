#' @title firstPackage
#'
#' @description This package computes the original Lord-Wingersky Algorithm for unidimensional IRT models, as well as the Lord-Wingersky Algorithm 2.0.
#'
#' @param thetaMin
#'
#' @param thetaMax
#'
#' @param nQuad
#'
#' @return popDist
#'
#' @examples normal_distribution(thetaMin=-5, thetaMax=5, nQuad=21)
#'
#' @export normal_distribution


normal_distribution <- function(nQuad,
                                thetaMin,
                                thetaMax) {
  qPoints <- seq(from = thetaMin,
                 to = thetaMax,
                 by = (thetaMax-thetaMin)/(nQuad-1))

  popDist <- exp(-(qPoints^2)/2)
  popDist <- popDist/sum(popDist)
  return(popDist)

}
