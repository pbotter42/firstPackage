#' @title firstPackage
#'
#' @description Computes a two dimensional trace surface
#'
#' @param theta_gen
#' theta values for the general dimension
#'
#' @param theta_spec
#' theta values for the specific dimension
#'
#' @param a_gen
#' descrimination parameter for the general dimension
#'
#' @param a_spec
#' descrimination parameter for the specific dimension(s)
#'
#' @param c
#' item threshold paramater
#'
#' @return ts
#'
#' @examples
#' mtr_example <- mtr(n_quad = 5,
#'                   theta_min = -2,
#'                   theta_max = 2,
#'                   a_gen = c(1.2,1.2,1,1,.8,.8),
#'                   a_spec = c(1,1,.8,.8,1.2,1.2),
#'                   c = c(-1,-.6,-.2,.2,.6,1))
#'
#' @export mar_tline

mtr <- function(n_quad,
                    theta_min,
                    theta_max,
                    a_gen,
                    a_spec,
                    c) {

  theta_gen <- quad_gen(n_quad, theta_min, theta_max)
  theta_spec <- quad_gen(n_quad, theta_min, theta_max)
  dist_2d <- norm_dist_2d(theta_gen, theta_spec)
  marg_2d <- marg_dist_2d(dist_2d)

  n_items <- length(a_gen)
  ts <- list()
  mar_tline <- list()
  for (i in 1:n_items) {
    ts[[i]] <- matrix(nrow = length(theta_spec),
                      ncol = length(theta_gen))
    mar_tline[[i]] <- rep(0, length(theta_gen))
    for (j in 1:length(theta_gen)) {
      for (k in 1:length(theta_spec)){
        ts[[i]][j,k] <- 1 / (1 + exp(-(a_gen[i]*theta_gen[j] + a_spec[i]*theta_spec[k] + (c[i]))))
        mar_tline[[i]][j] <- mar_tline[[i]][j]+ts[[i]][j,k]*dist_2d[j,k]
      }
      mar_tline[[i]][j] <- mar_tline[[i]][j]/marg_2d[j]
    }
  }
  return(mar_tline)
}
