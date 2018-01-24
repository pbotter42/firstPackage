#' @title firstPackage
#'
#' @description This package computes the original Lord-Wingersky Algorithm for unidimensional IRT models, as well as the Lord-Wingersky Algorithm 2.0.
#'
#' @param theta_gen
#'
#' @param theta_spec
#'
#' @return dist_m
#'
#' @examples normal_distribution(theta_gen = c(-2,-1,0,1,2),
#'  theta_spec = c(-2,-1,0,1,2))
#'
#' @export normal_distribution

norm_dist_2d <- function(theta_gen, theta_spec) {
  dist_m <- matrix(0,length(theta_gen),length(theta_gen))
  phi <- 0#toDo
  for (i in 1:length(theta_gen)) {
    for (j in 1:length(theta_spec)) {
      dist_m[i,j] <- exp(-0.5*(theta_gen[i]^2+theta_spec[j]^2-2*phi*theta_gen[i]*theta_spec[j])/(1-phi*phi))
    }
  }
  dist_m <- dist_m/sum(dist_m)
  return(dist_m)
}
