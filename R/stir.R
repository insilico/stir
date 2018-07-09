# package creation notes
# 
#library(devtools)
#install.packages("roxygen2")
#library(roxygen2)
#setwd()
#create("stir")
#document() # creates roxygen documentation
#install("stir")  # local
#devtools::use_vignette("STIRvignette")  # first time you create a vignette
# commit and push to github
#install_github("insilico/stir") # github
#library(stir)

#=========================================================================#
#' check.packages
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
check.packages <- function(pkg){
  # check.packages function: install and load multiple R packages.
  # Check to see if packages are installed. Install them if they are not, 
  # then load them into the R session.
  # https://gist.github.com/smithdanielle/9913897
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#=========================================================================#
#' attr.range
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
attr.range <- function(my.mat) {
  # compute denominator of the diff formula
  # for each attribute x (column) in my.mat, max(x) -  min(x)
  apply(as.matrix(my.mat), 2, function(x) {max(x) - min(x)})
}


#=========================================================================#
#' hamming.binary
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
hamming.binary <- function(X) {
  # compute hamming distance for a binary matrix
  D <- t(1 - X) %*% X
  D + t(D)
}


#=========================================================================#
#' sort.scores
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
sort.scores <- function(scores.vec){
  # sort attributes based on scores, important attributes on top
  sort(scores.vec, decreasing = TRUE)
}

#=========================================================================#
#' sort.pvalue
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
sort.pvalue <- function(pvalue.vec){
  # sort attributes based on pvalues, important attributes on top
  sort(pvalue.vec, decreasing = FALSE)
}

#=========================================================================#
#' diff.func
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
diff.func <- function(a, b, norm.fac = 1, metric){
  # compute the difference between two vectors elementwise
  if (metric %in% c("manhattan", "hamming")) val <- abs(a - b)/norm.fac
  if (metric == "euclidean") val <- abs(a - b)^2/norm.fac
  # more metrics here?
  val
}


#=========================================================================#
#' get.distance
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
get.distance <- function(attr.mat, metric){
  # Compute distance between two sample rows (instances)
  # based on all attributes, normalized by max-min
  # metric options: "euclidean", "maximum", "manhattan", "hamming"
  if (metric == "hamming"){
    distance.mat <- hamming.binary(attr.mat)
  } else {
    maxminVec <- attr.range(attr.mat)
    minVec <- apply(attr.mat, 2, function(x) {min(x)})
    attr.mat.centered <- t(attr.mat) - minVec
    attr.mat.scale <- t(attr.mat.centered / maxminVec)
    distance.mat <- as.matrix(dist(attr.mat.scale, method = metric))
  }
  distance.mat
}


#=========================================================================#
#' regular.ttest
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
regular.ttest.fn <- function(attr.idx, dat, class.idx = ncol(dat)){
  # perform t-test for each attribute with index attr.idx
  # assuming 2 classes
  dat[, class.idx] <- as.factor(dat[, class.idx])
  classes <- levels(dat[, class.idx])
  predictors.mat <- dat[, - class.idx]
  t.test(predictors.mat[dat[, class.idx] == classes[1], attr.idx], 
         predictors.mat[dat[, class.idx] == classes[2], attr.idx])$p.value
}


#=========================================================================#
#' find.neighbors
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
find.neighbors <- function(attr.mat, pheno.class, method, k, delta.vec = NULL){
  # create a matrix with num.samp rows, three columns
  # first column is sample Ri, second is Ri's hit, third is Ri's miss
  
  dist.mat <- get.distance(attr.mat, metric = metric)
  num.samp <- nrow(attr.mat)
  num.pair <- num.samp * (num.samp-1) / 2 # number of paired distances
  radius.surf <- sum(dist.mat)/(2*num.pair) # const r = mean(all distances)
  
  if (method == "relieff"){  
    neighbor.idx <- matrix(0, nrow = num.samp * k, ncol = 3)
    
    for (Ri in seq(1:num.samp)){ # for each sample Ri
      Ri.distances <- dist.mat[Ri,] # all distances to sample Ri
      Ri.nearest <- order(Ri.distances, decreasing = F) # closest to furthest
      Ri.nearest.mat <- cbind(Ri.nearest, pheno.class[Ri.nearest])
      Ri.nearest.hits <- Ri.nearest.mat[Ri.nearest.mat[,2] == Ri.nearest.mat[1,2],]
      Ri.nearest.misses <- Ri.nearest.mat[Ri.nearest.mat[,2] != Ri.nearest.mat[1,2],]
      Ri.kn.hits <- Ri.nearest.hits[2:(k+1), 1] # skip Ri self
      Ri.kn.misses <- Ri.nearest.misses[1:k, 1]
      
      # stack matrix of neighbor indices
      row.start <- (Ri-1)*k + 1
      row.end <- row.start + k - 1
      neighbor.idx[row.start:row.end, 1] <- rep(Ri, k)
      neighbor.idx[row.start:row.end, 2] <- Ri.kn.hits
      neighbor.idx[row.start:row.end, 3] <- Ri.kn.misses      
    }
    colnames(neighbor.idx) <- c("Ri_idx","hit_idx","miss_idx")
    hit.idx <- neighbor.idx[, c(1, 2)]
    miss.idx <- neighbor.idx[, c(1, 3)]
    
  } else {
    if (method == "surf") Ri.radius <- rep(radius.surf, num.samp) # use constance radius
    if (method == "multisurf"){
      if (is.null(delta.vec)) delta.vec <- sapply(1:num.samp, function(x) sd(dist.mat[-x, x]))
      Ri.radius <- colSums(dist.mat)/(num.samp - 1) - delta.vec/2 # use adaptive radius
    }
    
    hit.idx <- data.frame()
    miss.idx <- data.frame()
    
    for (Ri in seq(1:num.samp)){ # for each sample Ri
      Ri.distances <- sort(dist.mat[Ri,], decreasing = F)
      Ri.nearest <- Ri.distances[Ri.distances < Ri.radius[Ri]] # within the threshold
      # Ri.nearest.idx <- as.numeric(names(Ri.nearest))
      Ri.nearest.idx <- match(names(Ri.nearest), row.names(attr.mat))
      Ri.nearest.mat <- data.frame(Ri.nearest.idx, pheno.class[Ri.nearest.idx])
      Ri.nearest.hits <- Ri.nearest.mat[Ri.nearest.mat[,2] == Ri.nearest.mat[1,2],]
      Ri.nearest.misses <- Ri.nearest.mat[Ri.nearest.mat[,2] != Ri.nearest.mat[1,2],]
      if ((nrow(Ri.nearest.hits) > 1) & (nrow(Ri.nearest.misses) > 1)){
        hit.idx <- rbind(hit.idx, cbind(Ri, Ri.nearest.hits[-1, 1])) # exclude itself
        miss.idx <- rbind(miss.idx, cbind(Ri, Ri.nearest.misses[, 1]))
      } 
    }
    
    colnames(hit.idx) <- c("Ri_idx", "hit_idx")
    colnames(miss.idx) <- c("Ri_idx", "miss_idx")
  }
  list(hit.idx, miss.idx)
}


#=========================================================================#
#' make.factor
#'
#' Description
#'
#' @param test test
#' @return test test 
#' @examples
#' Example
#' @export
# multiply with 1/k (ReliefF) or 1/#hits and 1/#misses for each Ri
make.factor <- function(x) rep(1/x, x)


#=========================================================================#
#' stir
#'
#' Compute stir statistics for attributes given nearest hit/miss matrix
#'
#' @param attr.mat m x p matrix of m instances and p attributes 
#' @param neighbor.idx nearest hit/miss matrix from \code{find.neighbors}
#' @param method neighborhood method [\code{"multisurf"} (k=0) or \code{"relieff"} (specify k)]
#' @param m optional number of instances
#' @param k number of nearest hits/misses for \code{"relieff"} method (k=0 for \code{"multisurf"})
#' @param metric for distance matrix between instances (\code{"manhattan"} or \code{"euclidean"})
#' @param transform transformation of distances (\code{"None"}, \code{"sqrt"} or \code{"neglog"})
#' @return rs.list list: relief.score.df, arf.tstats.ordered, arf.fstats.ordered
#'
#' @examples
#' #See vignette("STIRvignette")
#' RF.method = "multisurf"
#' metric <- "manhattan"
#' neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.class, k = 0, method = RF.method)
#' results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
#' t_sorted_multisurf <- results.list$`STIR-t`
#' t_sorted_multisurf$attribute <- rownames(t_sorted_multisurf)
#'
#' @export
stir <- function(attr.mat, neighbor.idx, method, m = nrow(attr.mat), 
                   k, metric, transform = "None"){
  # simple implementation of the relieff algorithm
  # returns a named vector of attribute scores
  
  
  num.attr <- ncol(attr.mat)
  n.samp <- nrow(attr.mat)
  vecW <- rep(0, num.attr)
  names(vecW) <- colnames(attr.mat)
  range.vec <- attr.range(attr.mat)
  one_over_range <- 1/range.vec   # max - min of all attribtes
  one_over_m <- 1/n.samp          # m in paper notation
  
  # Initialize Relief-F T-stat and Anova F-stat
  arf.tstats <- matrix(0, nrow = num.attr, ncol = 3)
  arf.tstats.samp <- matrix(0, nrow = num.attr, ncol = 2)
  arf.fstats <- matrix(0, nrow = num.attr, ncol = 2)
  stir.tstats <- matrix(0, nrow = num.attr, ncol = 2)
  colnames(arf.tstats) <- c("t.stat", "t.pval", "cohen.d")
  colnames(arf.tstats.samp) <- c("t.stat", "t.pval")
  colnames(arf.fstats) <- c("F.stat", "F.pval")
  colnames(stir.tstats) <- c("t.stat.stir", "t.pval.stir")
  rownames(arf.tstats) <- colnames(attr.mat)
  rownames(arf.tstats.samp) <- colnames(attr.mat)
  rownames(arf.fstats) <- colnames(attr.mat)
  rownames(stir.tstats) <- colnames(attr.mat)
  
  hit.df <- neighbor.idx[[1]]
  miss.df <- neighbor.idx[[2]]  
  n.misses <- nrow(miss.df)    # num miss pairs (m*k for ReleifF)
  n.hits <- nrow(hit.df)       # num hit  pairs (m*k for ReleifF)
  #avg.k <- (n.misses + n.hits)/(2*n.samp)  # for subsampling
  
  hit.idx <- hit.df[, "hit_idx"]
  miss.idx <- miss.df[, "miss_idx"]
  Ri.hit.idx <- hit.df[, "Ri_idx"]
  Ri.miss.idx <- miss.df[, "Ri_idx"]
  
  nhit.a <- as.numeric(table(Ri.hit.idx))
  nmiss.a <- as.numeric(table(Ri.miss.idx))
  one_over_k_hits.fac <- unlist(lapply(nhit.a, make.factor))   # k_H_i vector; 1/k in general
  one_over_k_miss.fac <- unlist(lapply(nmiss.a, make.factor)) # k_M_i vector; 1/k in general
  
  for (attr.idx in seq(1, num.attr)){

    attr.val <- attr.mat[, attr.idx]
    hit.neighbors <- attr.val[hit.idx]
    miss.neighbors <- attr.val[miss.idx]
    Ri.hit.val <- attr.val[Ri.hit.idx]
    Ri.miss.val <- attr.val[Ri.miss.idx]

    # differencing vectors between Ri/hits and Ri/misses (1st and 2nd; 1st and 3rd)
    attr.diff.hit <- diff.func(hit.neighbors, Ri.hit.val, metric = metric) * one_over_range[attr.idx] # this is a vector
    attr.diff.miss <- diff.func(miss.neighbors, Ri.miss.val, metric = metric) * one_over_range[attr.idx] # this is a vector
    
    
    #######################################
    # 1. T-test
    #######################################   
    # assume equal variance estimate of t statistic:
    # otherwise:     
    # s.pooled <- sqrt(s.misses^2/n.misses + s.hits^2/n.hits)
    # t.stat.man <- (mu.misses - mu.hits)/(s.pooled)
    if (transform == "None"){
      attr.diff.trans.miss <- attr.diff.miss
      attr.diff.trans.hit <- attr.diff.hit      
    } else if (transform == "sqrt"){
      attr.diff.trans.miss <- sqrt(attr.diff.miss)
      attr.diff.trans.hit <- sqrt(attr.diff.hit)
    } else if (transform == "neglog"){
      attr.diff.trans.miss <- -log(attr.diff.miss)
      attr.diff.trans.hit <- -log(attr.diff.hit)
    }

    attr.diff.trans.miss <- attr.diff.trans.miss * one_over_k_miss.fac
    attr.diff.trans.hit <- attr.diff.trans.hit * one_over_k_hits.fac
    attr.diff.trans.miss.rf <- attr.diff.trans.miss # save for original relief scores later
    attr.diff.trans.hit.rf <- attr.diff.trans.hit # save for original relief scores later
    
    mu.misses <- sum(attr.diff.trans.miss)/n.samp  
    variance.misses <- sum( one_over_k_miss.fac * (attr.diff.miss-mu.misses)^2 )/(n.samp)
    s.misses <- sqrt(variance.misses)
    
    mu.hits <- sum(attr.diff.trans.hit)/n.samp  #bam
    variance.hits <- sum( one_over_k_hits.fac * (attr.diff.hit-mu.hits)^2 )/(n.samp)
    s.hits <- sqrt(variance.hits)
    
    s.adj <- 0
    s.pooled <- sqrt(((n.misses -1)*variance.misses + (n.hits - 1)*variance.hits)/(n.hits + n.misses - 2))
    t.stat.man <- (mu.misses - mu.hits)/(s.pooled*sqrt(1/n.misses + 1/n.hits))

    cohen.d <- t.stat.man*sqrt(1/n.misses+ 1/n.hits)
    r.miss <- variance.misses/n.misses
    r.hit <- variance.hits/n.hits
    # df <- floor((r.miss + r.hit)^2/(r.miss^2/(n.misses - 1) + r.hit^2/(n.hits - 1))) # degrees of freedom
    df <- n.hits + n.misses - 2
    man.pval <- pt(q = t.stat.man, df = df, lower.tail = F)
    arf.tstats[attr.idx, ] <- c(t.stat.man, man.pval, cohen.d)
    
    # or, use already built-in t-test in R's stats package
    # t.builtin <- t.test(attr.diff.trans.miss, attr.diff.trans.hit, alternative = "greater")
    # arf.tstats.builtin[attr.idx, ] <- c(t.builtin$statistic, t.builtin$p.value)
    
    #######################################
    # 1b. T-test Redo
    #######################################  

    # my.hit.df <- data.frame(Ri.hit.idx, attr.diff.hit)
    # summarised.t.hit <- my.hit.df %>% group_by(Ri.hit.idx) %>% summarise(avg = mean(attr.diff.hit))
    # 
    # my.miss.df <- data.frame(Ri.miss.idx, attr.diff.miss)
    # summarised.t.miss <- my.miss.df %>% group_by(Ri.miss.idx) %>% summarise(avg = mean(attr.diff.miss))
    # # mean(summarised.t.miss$avg) - mean(summarised.t.hit$avg)
    # 
    # t.test.stir <- t.test(summarised.t.miss$avg, summarised.t.hit$avg, alternative = "greater")
    # stir.tstats[attr.idx, ] <- c(t.test.stir$statistic, t.test.stir$p.value)
    # # mean(summarised.t.miss$avg) - mean(summarised.t.hit$avg)
      
      
      
      
    #######################################
    # 2. F-test
    #######################################   
    diff.aov.df <- data.frame(diffs = c(attr.diff.hit, attr.diff.miss),
                              hitmiss = c(rep("hit", n.hits), rep("miss", n.misses))) 
    aov.result <- summary(aov(diffs ~ hitmiss, data = diff.aov.df))
    arf.fstats[attr.idx,] <- as.numeric(aov.result[[1]][1, c("F value", "Pr(>F)")])
    
    #######################################
    # 3. ReliefF score
    #######################################   
    

    vecW[attr.idx] <- (sum(attr.diff.trans.miss.rf) - sum(attr.diff.trans.hit.rf)) 
    
  }
  

  relief.score <- vecW * one_over_m
  relief.score.df <- data.frame(attribute = names(relief.score), relief.score = relief.score)
  relief.score.ordered <- relief.score.df[order(relief.score.df[, "relief.score"], decreasing = F), ]
  
  # order results based on p-value
  arf.tstats.ordered <- data.frame(arf.tstats[order(arf.tstats[, "t.pval"], decreasing = F), ])
  arf.fstats.ordered <- data.frame(arf.fstats[order(arf.fstats[, "F.pval"], decreasing = F), ])
  # stir.tstats.ordered <- data.frame(stir.tstats[order(stir.tstats[, "t.pval.stir"], decreasing = F), ])
  
  # adjust p-values using Benjamini-Hochberg (default)
  arf.tstats.ordered$t.pval.adj <- p.adjust(arf.tstats.ordered[, "t.pval"])
  arf.fstats.ordered$F.pval.adj <- p.adjust(arf.fstats.ordered[, "F.pval"])
  # arf.tstats.ordered.samp$t.pval.adj <- p.adjust(arf.tstats.ordered.samp[, "t.pval"])
  # stir.tstats.ordered$t.pval.stir.adj <- p.adjust(stir.tstats.ordered[, "t.pval.stir"])
  rs.list <- list(relief.score.df, arf.tstats.ordered, arf.fstats.ordered)
  names(rs.list) <- c("OriRelief", "STIR-t", "STIR-F")
  return(rs.list)
}


#=========================================================================#
#' convert.pec.sim.to.inbix
#'
#' Reads a pEC simulated csv file, converts to inbix format and writes to working directory
#'
#' @param pEC.inputFile pEC-simulated data csv file name with path if needed 
#' @param inbix.file.prefix prefix of inbix filename for .num and .pheno files
#' @return none
#'
#' @examples
#' # setwd(...) # for output
#' infile <- "ARF_compare_1_multisurf_0.8_bias_No_k_0.1_pct.signals_1000_num.attr_100_num.samp.csv"
#' convert.pec.sim.to.inbix(infile,"simulated1")
#' # writes simulated1.num and simualted1.pheno
#'
#' @export
convert.pec.sim.to.inbix <- function(pEC.inputFile,inbix.file.prefix){
  
  ### write pEC simulated csv format to .num and .pheno format for inbix
  
  simdat <- read.csv(file=pEC.inputFile)
  # first column is X: subject id's case1, case2,...
  # class columns is class/status: 1, 1, -1, -1, ...
  nc <- ncol(simdat)
  num.attr <- nc - 2
  subjIDs <- as.character(simdat[,1])
  pheno <- as.numeric(simdat[, nc])  # class is -1/1, change to 1/2
  pheno[pheno==1]<-2
  pheno[pheno==-1]<-1
  predictors.mat <- simdat[,-nc]        # drop class (last) col
  predictors.mat <- predictors.mat[,-1]  # drop subj id (first) col
  attr.names <- colnames(predictors.mat)
  
  # for inbix .pheno. note: .pheno does not have header
  # write inbix .pheno
  phenosTable <- cbind(subjIDs, subjIDs, pheno)
  datasimInbixPhenoFile <- paste(inbix.file.prefix, ".pheno", sep="")
  write.table(phenosTable, datasimInbixPhenoFile, quote=F, sep="\t", 
              col.names=F, row.names=F)
  
  # .num has a header
  # write inbix numeric (.num) file/data set
  dataTable <- cbind(subjIDs, subjIDs, predictors.mat)
  colnames(dataTable) <- c("FID", "IID", attr.names)
  datasimInbixNumFile <- paste(inbix.file.prefix, ".num", sep="")
  write.table(dataTable, datasimInbixNumFile, quote=F, sep="\t", 
              col.names=T, row.names=F)
}
