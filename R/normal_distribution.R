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


normal_distribution <- function(n_quad,
                                theta_min,
                                theta_max) {
  q_points <- seq(from = theta_min,
                 to = theta_max,
                 by = (theta_max-theta_min)/(n_quad-1))

  dist_m <- exp(-(q_points^2)/2)
  dist_m <- dist_m/sum(dist_m)
  return(dist_m)

}


