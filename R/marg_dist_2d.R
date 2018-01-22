#' @title firstPackage
#'
#' @description This package computes the original Lord-Wingersky Algorithm for unidimensional IRT models, as well as the Lord-Wingersky Algorithm 2.0.
#'
#' @param n_quad
#'
#' @param theta_max
#'
#' @param theta_min
#'
#' @return q_points
#'
#' @examples quad_gen(n_quad=21,theta_min=-5, theta_max=5)
#'
#' @export marg_dist_2d
#'
marg_dist_2d <- function(dist_2d) {
  n_quad <- nrow(dist_2d)
  dist_m <- rep(0, n_quad)
  for (i in 1:n_quad) {
    for (j in 1:n_quad){
      dist_m[i] <- dist_m[i] + dist_2d[i,j]
    }
  }
  return(dist_m)
}

foo <- norm_dist_2d(21,-5, 5)
marg_dist_2d(foo)

