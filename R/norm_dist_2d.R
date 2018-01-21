#' @title firstPackage
#'
#' @description This package computes the original Lord-Wingersky Algorithm for unidimensional IRT models, as well as the Lord-Wingersky Algorithm 2.0.
#'
#' @param theta_min
#'
#' @param theta_max
#'
#' @param n_quad
#'
#' @return dist_m
#'
#' @examples normal_distribution(n_quad=21,theta_min=-5, theta_max=5)
#'
#' @export normal_distribution

norm_dist_2d <- function(n_quad,theta_min,theta_max) {
  dist_m <- matrix(0,n_quad,n_quad)
  theta_gen <- seq(from = theta_min,
                  to = theta_max,
                  by = (theta_max-theta_min)/(n_quad-1))# quad points for the general factor
  theta_spec <- seq(from = theta_min,
                 to = theta_max,
                 by = (theta_max-theta_min)/(n_quad-1))# quad points for the general factor
  phi <- 0#toDo
  for (i in 1:length(theta_gen)) {
    for (j in 1:length(theta_spec)) {
      dist_m[i,j] <- exp(-0.5*(theta_gen[i]^2+theta_spec[j]^2-2*phi*theta_gen[i]*theta_spec[j])/(1-phi*phi))
    }
  }
  dist_m <- dist_m/sum(dist_m)
  return(dist_m)
}
