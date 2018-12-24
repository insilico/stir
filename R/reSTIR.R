#=========================================================================#
#' nearestNeighbors
#'
#' Find nearest neighbors of each instance using relief.method
#' Used for regression stir (reSTIR) (no hits or misses specified in function).
#'
#' @param attr.mat m x p matrix of m instances and p attributes 
#' @param metric for distance matrix between instances (\code{"manhattan"} or \code{"euclidean"})
#' @param nbd.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]
#' @param k number of constant nearest hits/misses for \code{"relieff"}
#' @param sd.frac multiplier of the standard deviation of the distances when subtracting from average for SURF or multiSURF.
#' The multiSURF default is sd.frac=0.5: mean - sd/2 
#' @return  Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
#'
#' @examples
#' #See vignette("STIRvignette")
#' RF.method = "multisurf"
#' metric <- "manhattan"
#' neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.class, k = 0, method = RF.method)
#' results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
#' t_sorted_multisurf <- results.list$STIR_T
#' t_sorted_multisurf$attribute <- rownames(t_sorted_multisurf)
#'
#' @export
nearestNeighbors <- function(attr.mat, metric = "manhattan", nbd.method="multisurf", k=0, sd.vec = NULL, sd.frac = 0.5){
  # create a matrix with num.samp rows, two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  
  dist.mat <- get.distance(attr.mat, metric = metric)
  num.samp <- nrow(attr.mat)
  num.pair <- num.samp * (num.samp-1) / 2 # number of paired distances
  radius.surf <- sum(dist.mat)/(2*num.pair) # const r = mean(all distances)
  
  if (nbd.method == "relieff"){  
    Ri_NN.idxmat <- matrix(0, nrow = num.samp * k, ncol = 2)
    colnames(Ri_NN.idxmat) <- c("Ri_idx","NN_idx")
    for (Ri in seq(1:num.samp)){ # for each sample Ri
      Ri.distances <- dist.mat[Ri,] # all distances to sample Ri
      Ri.nearest <- order(Ri.distances, decreasing = F) # closest to furthest
      ## bam_add
      Ri.nearest.idx <- Ri.nearest[2:(k+1)] # skip Ri self
      # stack matrix of neighbor indices
      row.start <- (Ri-1)*k + 1
      row.end <- row.start + k - 1
      Ri_NN.idxmat[row.start:row.end, 1] <- rep(Ri, k)     # col of Ri's
      Ri_NN.idxmat[row.start:row.end, 2] <- Ri.nearest.idx  # col of knn's of Ri's
    }
  } else {
    if (nbd.method == "surf"){
      sd.const <- sd(dist.mat[upper.tri(dist.mat)])  
      # bam: orignal surf does not subtract sd-frac but should for fair multisurf comparison
      Ri.radius <- rep(radius.surf - sd.frac*sd.const, num.samp) 
    }
    if (nbd.method == "multisurf"){
      if (is.null(sd.vec)) sd.vec <- sapply(1:num.samp, function(x) sd(dist.mat[-x, x]))
      Ri.radius <- colSums(dist.mat)/(num.samp - 1) - sd.frac*sd.vec # use adaptive radius
    }
    
    # put each Ri's nbd in a list then rbind them at the end with do.call(rbind, List)
    # initialize list:
    Ri.nearestPairs.list <- vector("list",num.samp)
    for (Ri in seq(1:num.samp)){ # for each sample Ri
      Ri.distances <- sort(dist.mat[Ri,], decreasing = F)
      Ri.nearest <- Ri.distances[Ri.distances < Ri.radius[Ri]] # within the threshold
      Ri.nearest <- Ri.nearest[-1] # skip Ri self
      Ri.nearest.idx <- match(names(Ri.nearest), row.names(attr.mat))
      if (length(Ri.nearest.idx) > 1){ # if neighborhood not empty
        # cbind automatically repeats Ri
        Ri.nearestPairs.list[[Ri]] <- cbind(Ri, Ri.nearest.idx) 
      } 
    } # end for, now stack lists into matrix, do.call rbind
    
    Ri_NN.idxmat <- do.call(rbind, Ri.nearestPairs.list)
    colnames(Ri_NN.idxmat) <- c("Ri_idx","NN_idx")
  }
  # matrix of Ri's (first column) and their NN's (second column)
  return(Ri_NN.idxmat)
}
