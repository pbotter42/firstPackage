lw_2 <- function(n_quad, theta_min, theta_max,
                 a_gen,
                 a_spec,
                 c,
                 ic_index,
                 nr) {

  theta_gen <- quad_gen(n_quad, theta_min, theta_max)
  theta_spec <- quad_gen(n_quad, theta_min, theta_max)
  dist_2d <- norm_dist_2d(theta_gen, theta_spec)
  marg_2d <- marg_dist_2d(dist_2d)
  tlines <- comp_ts(theta_gen, theta_spec, a_gen, a_spec,
                    c, ic_index, nr)

  lw_lik <- list() # IC-lw step-item (1- lik, 2 - ss)
  lw_final_step <- list()
  for(k in 1:length(tlines)) { # i.e. for each item cluster

    lw_lik[[k]] <- list() # list at the level of item clusters
    lw_lik[[k]][[1]] <- list() # list for IC k lik cal iter 1
    lw_lik[[k]][[1]][[1]] <- list() # list for IC k lik cal iter 1 score 1
    lw_lik[[k]][[1]][[2]] <- list()# list for IC k lik cal iter 1 score 2
    lw_lik[[k]][[1]][[1]][[1]] <- matrix(unlist(tlines[[k]][1]),
                                         nrow = n_quad,
                                         ncol = n_quad)
    lw_lik[[k]][[1]][[1]][[2]] <- numeric()
    lw_lik[[k]][[1]][[1]][[2]] <- 1 # correct for p
    lw_lik[[k]][[1]][[2]][[1]] <- 1-matrix(unlist(tlines[[k]][1]),
                                           nrow = n_quad,
                                           ncol = n_quad)
    lw_lik[[k]][[1]][[2]][[2]] <- numeric()
    lw_lik[[k]][[1]][[2]][[2]] <- 0 # correct for p

    for(i in 2:length(tlines[[k]])) {
      lw_lik[[k]][[i]] <- c(lw_lik[[k]][[i-1]],
                            lw_lik[[k]][[i-1]])
      ts_num <- length(lw_lik[[k]][[i]])
      for(j in 1:ts_num) {
        if(j <= ts_num/2) {
          lw_lik[[k]][[i]][[j]][[1]] <- lw_lik[[k]][[i]][[j]][[1]]*tlines[[k]][[i]]
          lw_lik[[k]][[i]][[j]][[2]] <- lw_lik[[k]][[i]][[j]][[2]]+1
        }
        else{
          lw_lik[[k]][[i]][[j]][[1]] <- lw_lik[[k]][[i]][[j]][[1]]*(1-tlines[[k]][[i]])
        }
      }
    }
    max_likCal_iter <- length(lw_lik[[k]]) # keeping only the last iteration of the item cluster lw alg
    lw_lik[[k]] <- lw_lik[[k]][[max_likCal_iter]]
  }
  unique_SS_lik <- list() # list to contain lik for unique SS
  temp <- list() # The purpose of temp is to create an index so that sum scores can be looped through, bc the SS = 0 index cannot be refrenced in loops.

  # looping over lw_lik to get temp[[k]][[2]] to contain unique sum scores
  for(k in 1:length(lw_lik)) {
    temp[[k]] <- list()
    temp[[k]][[1]] <- numeric() # Vector for all Sum Scores across item clusters
    for(i in 1:length(lw_lik[[k]])) {
      temp[[k]][[1]][i] <- lw_lik[[k]][[i]][[2]]
    }
    temp[[k]][[2]] <- numeric()
    temp[[k]][[3]] <- numeric()
    temp[[k]][[2]] <- unique(temp[[k]][[1]])
    temp[[k]][[2]] <- sort(unique(temp[[k]][[1]]))
    temp[[k]][[3]] <- data.frame(SS = temp[[k]][[2]],
                                 index = seq(1, length(temp[[k]][[2]])))

  }

  # first simply create empty matrices(i.e. 0's) so that they can be used to sum lik later
  for(k in 1:length(temp[[k]])){
    unique_SS_lik[[k]] <- list()
    for(i in 1:length(temp[[k]][[2]])) { # indexing on SS +1 bc r cannot index on 0
      unique_SS_lik[[k]][[i]] <- list()
      unique_SS_lik[[k]][[i]][[1]] <- matrix(rep(0,
                                                 n_quad*n_quad),
                                             nrow = n_quad,
                                             ncol = n_quad)
      unique_SS_lik[[k]][[i]][[2]] <- numeric()
      unique_SS_lik[[k]][[i]][[2]] <- temp[[k]][[2]][i]

    }
  }

  # the loop below combines lik
  for(k in 1:length(lw_lik)) {
    for(i in 1:length(lw_lik[[k]])) {
      for(j in 1:nrow(temp[[k]][[3]])) {
        if(lw_lik[[k]][[i]][[2]] == temp[[k]][[3]][[j,"SS"]])
        {
          unique_SS_lik[[k]][[j]][[1]] <- unique_SS_lik[[k]][[j]][[1]]+lw_lik[[k]][[i]][[1]]
        }
      }
    }
  }

  for(k in 1:length(unique_SS_lik)) {
    for(i in 1:length(unique_SS_lik[[k]])) {
      unique_SS_lik[[k]][[i]][["mar_tLine"]] <- rep(0,n_quad) # creating a numeric vector that will be the marginal trace line
      for(QP1 in 1:n_quad) {
        for(QP2 in 1:n_quad) {
          unique_SS_lik[[k]][[i]][["mar_tLine"]][QP1] <- unique_SS_lik[[k]][[i]][["mar_tLine"]][QP1]+unique_SS_lik[[k]][[i]][[1]][QP1,QP2]*normal2d[QP1,QP2]
        }
        unique_SS_lik[[k]][[i]][["mar_tLine"]][QP1] <-  unique_SS_lik[[k]][[i]][["mar_tLine"]][QP1]/marg_2d[QP1]
      }
    }
  }

  unique_SS_lik2 <- list()
  for(k in 1:length(unique_SS_lik)) {
    unique_SS_lik2[[k]] <- list()
    for(i in 1:length(unique_SS_lik[[k]])) {
      unique_SS_lik2[[k]][[i]] <- list()
      unique_SS_lik2[[k]][[i]][["mar_tLine"]] <- unique_SS_lik[[k]][[i]]$mar_tLine
      unique_SS_lik2[[k]][[i]][["ic_SS"]] <- unique_SS_lik[[k]][[i]][[2]]
    }
  }
  lw_iter <- list()
  lw_iter[[1]] <- unique_SS_lik2[[1]]# the list for iteration1 of the part 2 recursions
  for(k in 2:length(unique_SS_lik2)) { # this loops through the 2nd item cluster to the last
    lw_iter[[k]] <- rep(lw_iter[[k-1]], length(unique_SS_lik2[[k]])) # notice how lw_iter is being repeated by the lengh of k
    n_iter <- 0
    for(i in 1:length(unique_SS_lik2[[k]])) {
      for(j in 1:length(lw_iter[[k-1]])) {
        n_iter <- n_iter+1
        lw_iter[[k]][[n_iter]]$mar_tLine <- lw_iter[[k]][[n_iter]]$mar_tLine*unique_SS_lik2[[k]][[i]]$mar_tLine # summing the lik
        lw_iter[[k]][[n_iter]]$ic_SS <- lw_iter[[k]][[n_iter]]$ic_SS+unique_SS_lik2[[k]][[i]]$ic_SS # summing the sum scores
      }
    }
  }
  lw_final_iter <- lw_iter[[length(lw_iter)]] # Returning only the final iteration of combining lik within item cluster

  lw2.0_Lik <- matrix(nrow = length(lw_final_iter),
                      ncol = n_quad)# Matrix for all the lik before they are combined via SS
  lw2.0_SS <- numeric() # vector to contain the SS

  for(k in 1:length(lw_final_iter)) {
    lw2.0_Lik[k,] <- lw_final_iter[[k]]$mar_tLine
    lw2.0_SS[k] <- lw_final_iter[[k]]$ic_SS
  }

  qPoint_names <- vector() # vector to store quad points names
  for(i in 1:n_quad) {
    qPoint_names[i] <- paste("Quad point", i, sep = " ")
  }
  colnames(lw2.0_Lik) <- qPoint_names
  lw2.0_Lik_SS <- cbind(lw2.0_Lik, lw2.0_SS) # contains the lik and SS
  lw2.0_Lik_SS <- as.data.frame(lw2.0_Lik_SS)
  aggregated_SS <- aggregate(.~lw2.0_SS, lw2.0_Lik_SS, sum)

  return(aggregated_SS)
}

Lis_example <- lw_2(n_quad = 5,
                    theta_min = -2,
                    theta_max = 2,
                    a_gen = c(1.2,1.2,1,1,.8,.8),
                    a_spec = c(1,1,.8,.8,1.2,1.2),
                    c = c(-1,-.6,-.2,.2,.6,1),
                    ic_index = c(1,1,2,2,3,3),
                    nr = 2)

Lis_example-f2
Lis_example2 <- lw_2(n_quad = 5,
                    theta_min = -2,
                    theta_max = 2,
                    a_gen = c(1.2,1.2,1,1,.8,.8),
                    a_spec = c(1,1,.8,.8,1.2,1.2),
                    c = c(-1,-.6,-.2,.2,.6,1),
                    ic_index = c(1,1,2,2,3,3),
                    nr = 2)

plot(c(-2,-1,0,1,2), Lis_example[1,2:6])
lines(c(-2,-1,0,1,2), Lis_example[2,2:6])
lines(c(-2,-1,0,1,2), Lis_example[3,2:6])
lines(c(-2,-1,0,1,2), Lis_example[4,2:6])
lines(c(-2,-1,0,1,2), Lis_example[5,2:6])
lines(c(-2,-1,0,1,2), Lis_example[6,2:6])
lines(c(-2,-1,0,1,2), Lis_example[7,2:6])

