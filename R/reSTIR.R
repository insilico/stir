#=========================================================================#
#' stirDiff
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
stirDiff <- function(a, b, type = "manhattan", norm.fac = 1){
  # compute the difference between two vectors elementwise
  if (type=="euclidean"){ # numeric
    val <- abs(a - b)^2/norm.fac
  } else if (type=="allele-sharing"){ # snps
    val <- abs(a-b)/2
  } else{ # manhattan, numeric
    val <- abs(a - b)/norm.fac
  }
  return(val)
}

#=========================================================================#
#' stirDistances
#'
#' Note: Should probalby standardize data before manhattan and euclidean?
#'
#'
#' @param attr.mat m x p matrix of m instances and p attributes 
#' @param metric for distance matrix between instances (default: \code{"manhattan"}, \code{"euclidean"}, 
#' \code{"relief-scaled-manhattan"}, \code{"relief-scaled-euclidean"}, \code{"allele-sharing-manhattan"}).
#' @examples
#' Example
#' @export
stirDistances <- function(attr.mat, metric="manhattan"){
  # Compute distance matrix between all samples (rows)
  # reSTIR default is numeric manhattan ("manhattan"), max-min scaling is not needed for stir
  if (metric == "hamming"){
    distance.mat <- hamming.binary(attr.mat)
  } else if (metric == "allele-sharing-manhattan"){
    # allele-sharing-manhattan, AM for SNPs
    attr.mat.scale <- attr.mat / 2
    distance.mat <- as.matrix(dist(attr.mat.scale, method = "manhattan"))
  } else if (metric == "relief-scaled-manhattan"){
    # value of metric, euclidean, manhattan or maximum
    maxminVec <- attr.range(attr.mat)
    minVec <- apply(attr.mat, 2, function(x) {min(x)})
    attr.mat.centered <- t(attr.mat) - minVec
    attr.mat.scale <- t(attr.mat.centered / maxminVec)
    distance.mat <- as.matrix(dist(attr.mat.scale, method = "manhattan"))
  } else if (metric == "relief-scaled-euclidean"){
    # value of metric, euclidean, manhattan or maximum
    maxminVec <- attr.range(attr.mat)
    minVec <- apply(attr.mat, 2, function(x) {min(x)})
    attr.mat.centered <- t(attr.mat) - minVec
    attr.mat.scale <- t(attr.mat.centered / maxminVec)
    distance.mat <- as.matrix(dist(attr.mat.scale, method = "euclidean"))
  } else if (metric=="euclidean"){
    distance.mat <- as.matrix(dist(attr.mat, method = "euclidean"))
  } else {
    distance.mat <- as.matrix(dist(attr.mat, method = "manhattan"))
  }
  distance.mat
}

#=========================================================================#
#' nearestNeighbors
#'
#' Find nearest neighbors of each instance using relief.method
#' Used for regression stir (reSTIR) (no hits or misses specified in function).
#' Also used in reSTIR for case/control, but hit/miss is used in reSTIR function. 
#'
#' @param attr.mat m x p matrix of m instances and p attributes 
#' @param metric used in stirDistances for distance matrix between instances, default: \code{"manhattan"} (numeric)
#' @param nbd.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]
#' @param k number of constant nearest hits/misses for \code{"relieff"}. 
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) 
#' @param sd.frac multiplier of the standard deviation of the distances when subtracting from average for SURF or multiSURF.
#' The multiSURF default is sd.frac=0.5: mean - sd/2 
#' @return  Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
#'
#' @examples
#' 
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, metric="manhattan", nbd.method = "multisurf", sd.frac = 0.5)
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, metric="manhattan", nbd.method = "relieff")
#'
#' @export
nearestNeighbors <- function(attr.mat, metric = "manhattan", nbd.method="multisurf", sd.vec = NULL, sd.frac = 0.5, k=0){
  # create a matrix with num.samp rows, two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  
  dist.mat <- stirDistances(attr.mat, metric = metric)
  num.samp <- nrow(attr.mat)
  
  if (nbd.method == "relieff"){  
    if (k==0){ # if no k specified or value 0
    # replace k with the theoretical expected value for SURF (close to multiSURF)
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    k <- floor((num.samp-1)*(1-erf(sd.frac/sqrt(2)))/2)
    }
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
      num.pair <- num.samp * (num.samp-1) / 2 # number of paired distances
      radius.surf <- sum(dist.mat)/(2*num.pair) # const r = mean(all distances)
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

#=========================================================================#
#' diffRegression
#'
#' Utility to run regression for a phenotype diff vector and one attribute diff vector. 
#' Organize regression statistics into a vector.
#'
#' @param pheno.diffs outcome variable as vector of diffs 
#' @param predictor.diffs one predictor varialbe as vector of diffs 
#' @param regression.type (\code{"linear"}) todo: \code{"logistic"} and covariates
#' @return vector or regression stats to put into list for reSTIR and combine into matrix
#'
#' @examples
#'
#' @export
diffRegression <- function(pheno.diffs, predictor.diffs, regression.type="lm") {
  # 
  if (regression.type=="lm"){
  fit <- summary(lm(pheno.diffs ~ predictor.diffs))
  }
  stats.vec <- c(fit$coefficients[1,3], # beta_hat_0, intercept
              fit$coefficients[2,3], # beat_hat_a, standardize beta for attribute
              fit$coefficients[1,4], # p-value for intercept
              fit$coefficients[2,4], # p-value for attribute beta
              fit$r.squared,         # R^2 of fit
              fit$fstatistic[1],     # F-stat and next is its p-value
              1 - pf(fit$fstatistic[1], fit$fstatistic[2], fit$fstatistic[3])
  )
  return(stats.vec)
}

#=========================================================================#
#' reSTIR
#'
#' regression-based STIR statistics for quantitative or case-control phenotypes
#' 
#'
#' @param pheno.vec length-m numeric outcome vector for linear regression, factor for logistic regression 
#' @param attr.mat m x p matrix of m instances and p attributes, include attr names as colnames 
#' @param neighbor.pairs.idx nearest hit/miss matrices, output from \code{find.neighbors}
#' @param attr.diff.type diff type for attributes (\code{"manhattan"} or \code{"euclidean"} for numeric)
#' @param pheno.diff.type diff type for phenotype (\code{"manhattan"} or \code{"euclidean"} for numeric)
#' @param regression.type (\code{"linear"} or \code{"logistic"})
#' @return reSTIR.results.df: reSTIR regression coefficients and p-values for each attribute
#'
#' @examples
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, metric="manhattan", nbd.method = "multisurf", sd.frac = 0.5)
#' results.mat <- reSTIR(pheno.vec, predictors.mat, neighbor.pairs.idx, attr.diff.type="manhattan", pheno.diff.type="manhattan", regression.type="linear")
#'
#' @export
reSTIR <- function(pheno.vec, attr.mat, neighbor.pairs.idx, attr.diff.type="manhattan", pheno.diff.type="manhattan", regression.type="linear"){

  num.attr <- ncol(attr.mat)
  num.samp <- nrow(attr.mat)
  
  # run reSTIR, each attribute is a list, then we do.call rbind to a matrix
  reSTIR.stats.list <- vector("list",num.samp)
  for (attr.idx in seq(1, num.attr)){
    attr.vals <- attr.mat[, attr.idx]
    Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1]]
    NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2]]
    attr.diff.vec <- stirDiff(Ri.attr.vals, NN.attr.vals)
    
    Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
    NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
    pheno.diff.vec <- stirDiff(Ri.pheno.vals, NN.pheno.vals)
    
    # utility function
    reSTIR.stats.list[[attr.idx]] <- diffRegression(pheno.diff.vec, attr.diff.vec)
  }

  reSTIR.stats.attr_ordered.mat <- do.call(rbind, reSTIR.stats.list)
  # todo: adjusted p-value of attribute beta and sort
  colnames(reSTIR.stats.attr_ordered.mat) <- c("B0", "Ba", "pval.0", "pval.a", "R.sqr", "Fstat", "Fstat.pval")
  
  if (!is.null(colnames(test.mat))){
    # add attribute names to stats/results matrix if the data matrix contains them
    rownames(reSTIR.stats.attr_ordered.mat) <- colnames(predictors.mat)
  }
  
  # order by attribute p-value
  reSTIR.results.ordered.mat <- reSTIR.stats.attr_ordered.mat[order(reSTIR.stats.attr_ordered.mat[, "pval.a"], decreasing = F), ]
  reSTIR.stats.df <- data.frame(reSTIR.results.ordered.mat)
  
  # order results based on p-value
  #arf.tstats.ordered <- data.frame(arf.tstats[order(arf.tstats[, "t.pval"], decreasing = F), ])
  #arf.fstats.ordered <- data.frame(arf.fstats[order(arf.fstats[, "F.pval"], decreasing = F), ])
  # stir.tstats.ordered <- data.frame(stir.tstats[order(stir.tstats[, "t.pval.stir"], decreasing = F), ])
  
  # adjust p-values using Benjamini-Hochberg (default)
  #arf.tstats.ordered$t.pval.adj <- p.adjust(arf.tstats.ordered[, "t.pval"])
  #arf.fstats.ordered$F.pval.adj <- p.adjust(arf.fstats.ordered[, "F.pval"])
  # arf.tstats.ordered.samp$t.pval.adj <- p.adjust(arf.tstats.ordered.samp[, "t.pval"])
  # stir.tstats.ordered$t.pval.stir.adj <- p.adjust(stir.tstats.ordered[, "t.pval.stir"])
  #rs.list <- list(OriRelief=relief.score.df, STIR_T=arf.tstats.ordered, STIR_F=arf.fstats.ordered)
  #names(rs.list) <- c("OriRelief", "STIR-t", "STIR-F")
  return(reSTIR.results.df)
}
