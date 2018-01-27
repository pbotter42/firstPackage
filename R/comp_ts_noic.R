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
#' location parameter for the general dimension
#'
#' @param a_spec
#' location parameter for the specific dimension(s)
#'
#' @param c
#' Guessing paramater
#'
#' @param ic_index
#' Numeric vector that delineates the correspondence between item and
#' item cluster.
#'
#' @param nr
#' number of possible response options
#'
#' @return ts_ic
#'
#' @examples
#' example_ts <- comp_ts(theta_gen = c(-2,-1,0,1,2),
#' theta_spec = c(-2,-1,0,1,2),
#' a_gen = c(1.2,1.2,1,1,.8,.8),
#' a_spec = c(1,1,.8,.8,1.2,1.2),
#' c = c(-1,-.6,-.2,.2,.6,1),
#' ic_index = c(1,1,2,2,3,3),
#' nr = 2)
#'
#' @export comp_ts

comp_ts <- function(theta_gen,
                    theta_spec,
                    a_gen,
                    a_spec,
                    c,
                    ic_index,
                    nr) {
  n_items <- length(a_gen)
  ts <- list()
  for (i in 1:n_items) {
    ts[[i]] <- matrix(nrow = length(theta_spec),
                      ncol = length(theta_gen))
    for (j in 1:length(theta_gen)) {
      for (k in 1:length(theta_spec)){
        ts[[i]][j,k] <- 1 / (1 + exp(-(a_gen[i]*theta_gen[j] + a_spec[i]*theta_spec[k] + (c[i]))))
      }
    }
  }
  ts_ic <- list() # List similar to ts, only that the trace surfaces will be indexed by item clusters
  ic_unique <- unique(ic_index) # Unique item cluster values
  item_index <- rep(seq(1,nr, by = 1), length(unique(ic_index)))
  n_iter <- 0
  for(i in 1:length(ic_unique)) {
    ts_ic[[i]] <- list()
    TF_index <- ic_index == ic_unique[i]
    for(j in 1:length(TF_index[TF_index == TRUE])) {
      n_iter <- n_iter+1
      ts_ic[[i]][[j]] <- matrix(nrow = length(theta_spec),
                                ncol = length(theta_gen))
      ts_ic[[ic_unique[i]]][[j]] <- ts[[n_iter]] #p
    }
  }
  return(ts_ic)
}
