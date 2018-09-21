STIR (STatistical Inference Relief) Example
================
Brett McKinney and Trang Le
2018-07-18

Install STIR and privateEC:
---------------------------

``` r
knitr::opts_chunk$set(echo = TRUE, warning=FALSE)
rm(list = ls())

if (!("devtools" %in% installed.packages()[,"Package"])){
  install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}
library(devtools)

if (!("privateEC" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/privateEC") # build_vignettes = TRUE)
}
if (!("stir" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/stir", build_vignettes = TRUE)
}
library(privateEC)  # used to simulate data
library(stir)

# load other helper packages
packages <- c("ggplot2", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
check.packages(packages)  # helper function from STIR
```

Simulate data with privateEC
----------------------------

``` r
letsSimulate <- T   # F to use previously simulated data
class.lab <- "class"
writeData <- F  # usually the same as letsSimulate
writeResults <- F

num.samp <- 100
num.attr <- 1000
pct.signals <- 0.1
bias <- 0.4
#sim.type <- "mainEffect"
sim.type <- "interactionErdos"

pec_simFile <- paste("pec_simulated", sim.type, "bias", bias, 
                             "pct.signals", pct.signals,
                             "num.attr", num.attr, "num.samp", num.samp, sep = "_")
pec_simFile <- paste(pec_simFile,".csv",sep="")

if (letsSimulate == TRUE){
    sim.data <- createSimulation(num.samples = num.samp, num.variables = num.attr,
                                 pct.signals = pct.signals, pct.train = 1/2, pct.holdout = 1/2, 
                                 bias = bias, sim.type = sim.type, verbose = FALSE)
  dat <- rbind(sim.data$train, sim.data$holdout)
  predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
} else { # optional: use provided data
  dat <- read.csv(pec_simFile)
  dat <- dat[,-1] # written file has first X column with subject names
  predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
}

dat[, class.lab] <- as.factor(dat[, class.lab]) 
pheno.class <- dat[, class.lab]
attr.names <- colnames(predictors.mat)
num.samp <- nrow(dat)

if (writeData == TRUE){
  write.csv(dat, file = pec_simFile)
}
```

### Run STIR-multiSURF:

``` r
RF.method = "multisurf"
metric <- "manhattan"
# let k=0 because multisurf does not use k
neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.class, k = 0, method = RF.method)
results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
t_sorted_multisurf <- results.list$STIR_T[, -3]  # remove cohen-d
colnames(t_sorted_multisurf) <- paste(c("t.stat", "t.pval", "t.pval.adj"), "stir", sep=".")
(t_sorted_multisurf[1:10,])
```

    ##          t.stat.stir  t.pval.stir t.pval.adj.stir
    ## simvar53   10.405782 2.688475e-25    2.688475e-22
    ## simvar69    8.979117 2.179423e-19    2.177244e-16
    ## simvar19    8.681967 2.931270e-18    2.925408e-15
    ## simvar49    8.402737 3.128784e-17    3.119397e-14
    ## simvar10    8.226849 1.339647e-16    1.334288e-13
    ## simvar48    7.678817 1.034552e-14    1.029379e-11
    ## simvar61    7.661820 1.178561e-14    1.171490e-11
    ## simvar12    7.621981 1.597942e-14    1.586757e-11
    ## simvar91    7.615994 1.672525e-14    1.659145e-11
    ## simvar57    7.369054 1.066236e-13    1.056640e-10

``` r
t_sorted_multisurf$attribute <- rownames(t_sorted_multisurf) # adds a column for merge
```

### Run STIR-ReliefF constant *k* = ⌊(*m* − 1)/6⌋:

ReliefF with *k* = ⌊(*m* − 1)/6⌋ (where m is the number of samples) is similar to multiSURF:

``` r
t_sorted_relieff <- list()
i <- 0
RF.method = "relieff"
k <- floor(num.samp/6)  # k=m/6 should be similar to MultiSURF
i <- i+1  # if you want to use k for loop
neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.class, k = k, method = RF.method)
results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
t_sorted_relieff[[i]] <- results.list$STIR_T[, -3]
colnames(t_sorted_relieff[[i]]) <- paste(c("t.stat", "t.pval", "t.pval.adj"), k, sep=".")
(t_sorted_relieff[[i]][1:10,])
```

    ##          t.stat.16    t.pval.16 t.pval.adj.16
    ## simvar53 10.436254 2.132324e-25  2.132324e-22
    ## simvar19  8.747637 1.729603e-18  1.727874e-15
    ## simvar69  8.464474 1.931650e-17  1.927787e-14
    ## simvar49  7.929971 1.501494e-15  1.496989e-12
    ## simvar10  7.895881 1.964286e-15  1.956429e-12
    ## simvar61  7.738830 6.678804e-15  6.645410e-12
    ## simvar91  7.629877 1.540153e-14  1.530912e-11
    ## simvar57  7.573471 2.363406e-14  2.346862e-11
    ## simvar12  7.543726 2.958625e-14  2.934956e-11
    ## simvar1   7.515995 3.645166e-14  3.612360e-11

``` r
t_sorted_relieff[[i]]$attribute <- rownames(t_sorted_relieff[[i]])
t_sorted_relieff[[i+1]] <- t_sorted_multisurf
```

### Standard t-test:

``` r
sort.pvalue <- function(pvalue.vec){
  # sort attributes based on pvalues, important attributes on top
  sort(pvalue.vec, decreasing = FALSE)
}
regular.ttest.results <- sapply(1:ncol(predictors.mat), regular.ttest.fn, dat = dat)
names(regular.ttest.results) <- colnames(predictors.mat)
regular.ttest.sorted <- sort.pvalue(regular.ttest.results)
regular.t.padj <- data.frame(regT.padj = p.adjust(regular.ttest.sorted))
```

### Aggregate results of STIR with ReleifF and STIR with MultiSURF and regular t-test:

``` r
final.mat <- Reduce(function(x, y) merge(x, y, by = "attribute", sort = F), t_sorted_relieff)
#final.mat <- reshape::merge_all(t_sorted_relieff)

# Are the columns sorted separately? There is only one column of attribute names 
# View(final.mat[1:15,],"Resutls: First 15 Rows")  # View has a problem with Rmarkdown

if (writeResults == T){
write.csv(final.mat,file="final.mat.csv")
}
```

Plot STIR significance of attributes:
-------------------------------------

``` r
rownames(final.mat) <- final.mat$attribute
pval.df <- final.mat[attr.names, ]

pval.melt <- melt(pval.df[, c("attribute", "t.pval.adj.stir", "t.pval.adj.16")], id.vars = 1)
levels(pval.melt$variable) <- c("multiSURF", "ReliefF, k=16")
pval.melt$value <- -log(pval.melt$value, 10)
pval.melt$value[pval.melt$value >10] <- 10

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
t4 <- ggplot(pval.melt, aes(x = attribute, y = value, group = variable, color = variable)) + 
  ylim(c(-0.2,11))+
  geom_point(alpha = 0.7, position = position_jitter(w = 0, h = 0.2)) + 
  geom_vline(xintercept = 100, linetype = 2, color = "grey") + 
  labs(y = "-Log10(p-value)", x = "Features (in data order)", title = "Significance of attributes") + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), axis.text.x=element_blank()) + 
  scale_color_manual(values = cbPalette[2:3]) +
  geom_hline(yintercept = -log(0.05, 10), linetype = 4, color = "grey") 
```

Plot of -log10(p-values) of attributes. Attributes are in their original order from the data, but the significant attributes tend to be on the left because the simulated functional attributes were targeted to be first. Thus, attributes to the left of the vertical dashed line are targeted as *functional* or *predictive* in the simulation. However, for interactions, some attributes on the right may be functional due to network co-expression. (Note: p-values less than *e*<sup>−10</sup> are plotted as *e*<sup>−10</sup> for scaling. Points are slightly jittered vertically to show results of both methods.)

``` r
show(t4)
```

![](STIRexample_files/figure-markdown_github/unnamed-chunk-7-1.png)
