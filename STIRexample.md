---
title: "A STIR Example"
author: "Brett McKinney and Trang Le"
date: '2018-06-29'
output:
  html_document:
    keep_md: yes
    knitr:
      clean: no
      run_pandoc: no
---

## Set up: 

Load packages, set parameters:


```r
#knitr::opts_chunk$set(echo = TRUE)
rm(list = ls())
#source('utilFuncs.R')
library(devtools)
```

```
## Warning: package 'devtools' was built under R version 3.4.3
```

```r
if (!("privateEC" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/privateEC")
}
#if (!("stir" %in% installed.packages()[,"Package"])){
#  devtools::install_github("insilico/stir")
#}
install_github("insilico/stir")
```

```
## Downloading GitHub repo insilico/stir@master
## from URL https://api.github.com/repos/insilico/stir/zipball/master
```

```
## Installing stir
```

```
## '/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file  \
##   --no-environ --no-save --no-restore --quiet CMD INSTALL  \
##   '/private/var/folders/s7/k42hr_yn2kz082lv_lnr4t9d4y8h_3/T/Rtmpd1Rv8x/devtools133271a7167b5/insilico-STIR-125d860'  \
##   --library='/Users/brett-mckinney/Library/R/3.4/library'  \
##   --install-tests
```

```
## 
```

```r
library(stir)
#install_github("insilico/privateEC")  # needed for simulations

library(reshape2)
```

```
## Warning: package 'reshape2' was built under R version 3.4.3
```

```r
library(ggplot2)
# load necessary packages
packages <- c("ggplot2", "devtools", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
check.packages(packages)
```

```
## Loading required package: CORElearn
```

```
## Warning: package 'CORElearn' was built under R version 3.4.4
```

```
## Loading required package: dplyr
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```
## Loading required package: pROC
```

```
## Warning: package 'pROC' was built under R version 3.4.4
```

```
## Type 'citation("pROC")' for a citation.
```

```
## 
## Attaching package: 'pROC'
```

```
## The following objects are masked from 'package:stats':
## 
##     cov, smooth, var
```

```
## Loading required package: plotROC
```

```
## 
## Attaching package: 'plotROC'
```

```
## The following object is masked from 'package:pROC':
## 
##     ggroc
```

```
##   ggplot2  devtools CORElearn  reshape2     dplyr      pROC   plotROC 
##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE
```

```r
library(privateEC)  # simulate data
```

```
## Loading privateEC
```

```r
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

letsSimulate <- T
metric <- "manhattan"
RF.method <- "multisurf"
writeResults <- T
class.lab <- "class"
writeData <- T
```

## Simulate data:


```r
if (letsSimulate == TRUE){
  n.samp <- 100
  num.samp <- n.samp
  num.attr <- 1000
  pct.signals <- 0.1
  bias <- 0.4
  run.permute <- TRUE
  is.main <- F  # simulate main effects or interactions

  if (is.main){ # make this main effect sim if TRUE
    sim.data <- createSimulation(num.samples = n.samp, num.variables = num.attr,
                                 pct.signals = pct.signals, pct.train = 1/2, pct.holdout = 1/2, 
                                 bias = bias, sim.type = "mainEffect", verbose = FALSE)
  } else { # interaction simulation
    sim.data <- createSimulation(num.samples = n.samp, num.variables = num.attr,
                                 pct.signals = pct.signals, pct.train = 1 / 2, pct.holdout = 1 / 2,
                                 bias = bias, sim.type = "interactionErdos", verbose = FALSE)
  }
  dat <- rbind(sim.data$train, sim.data$holdout)
} else {
  dat <- read.csv("ARF_compare_1_multisurf_0.8_bias_No_k_0.1_pct.signals_1000_num.attr_100_num.samp.csv")
  n.samp <- nrow(dat)
}

dat[, class.lab] <- as.factor(dat[, class.lab]) 
pheno.class <- dat[, class.lab]
predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
attr.names <- colnames(predictors.mat)
num.samp <- nrow(dat)
neighbor.idx <- find.neighbors(predictors.mat, method = RF.method, pheno.class, k = k)

if (writeData == TRUE){
  write.csv(dat, file = paste("ARF_compare", RF.method, bias, "bias", 
                             pct.signals, "pct.signals",
                             num.attr, "num.attr", num.samp, "num.samp.csv", sep = "_"))
}
```

### Run multiSURF:


```r
RF.method = "multisurf"
metric <- "manhattan"
# bam: I'm letting k=0 because multisurf
neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.class, k = 0, method = RF.method)
results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
# t_observed_mat <- results.list[[4]]
t_sorted_multisurf <- results.list$`STIR-t`
# vecW_observed <- results.list[[1]]
t_sorted_multisurf$attribute <- rownames(t_sorted_multisurf)
```

### Run ReliefF:

Run an example ReliefF with $k=\lfloor(m-1)/6\rfloor$:


```r
t_sorted_relieff <- list()
i <- 0
RF.method = "relieff"
# for (k in c(5, 10, floor(n.samp/6), floor(n.samp/3), floor((n.samp-1)/2))){
k <- floor(n.samp/6)
i <- i+1
neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.class, k = k, method = RF.method)
results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
t_sorted_relieff[[i]] <- results.list$`STIR-t`[, -3]
colnames(t_sorted_relieff[[i]]) <- paste(c("t.stat", "t.pval", "t.pval.adj"), k, sep=".")
t_sorted_relieff[[i]]$attribute <- rownames(t_sorted_relieff[[i]])
t_sorted_relieff[[i+1]] <- t_sorted_multisurf
```

### Run standard t-test:


```r
regular.ttest.results <- sapply(1:ncol(predictors.mat), regular.ttest.fn, dat = dat)
names(regular.ttest.results) <- colnames(predictors.mat)
regular.ttest.sorted <- sort.pvalue(regular.ttest.results)
regular.t.padj <- data.frame(regT.padj = p.adjust(regular.ttest.sorted))
```

### Aggregate results:


```r
final.mat <- Reduce(function(x, y) merge(x, y, by = "attribute", sort = F), t_sorted_relieff)
# final.mat <- reshape::merge_all(t_sorted_relieff)
write.csv(final.mat,file="final.mat.csv")
```


## Plot results:

Plot the significance of attributes. p-values $< e^{-10}$ are plotted as $< e^{-10}$ for better visual scale. Attributes to the left of the vertical dash line are targeted as *functional* or *predictive* in the simulation. (Dots are slightly jittered vertically to show both methods' results.)


```r
rownames(final.mat) <- final.mat$attribute
pval.df <- final.mat[attr.names, ]

pval.melt <- melt(pval.df[, c("attribute", "t.pval.adj", "t.pval.adj.16")], id.vars = 1)
levels(pval.melt$variable) <- c("multiSURF", "ReliefF, k=16")
pval.melt$value <- -log(pval.melt$value, 10)
pval.melt$value[pval.melt$value >10] <- 10

t4 <- ggplot(pval.melt, aes(x = attribute, y = value, group = variable, color = variable)) + 
  ylim(c(-0.2,11))+
  geom_point(alpha = 0.7, position = position_jitter(w = 0, h = 0.2)) + 
  geom_vline(xintercept = 100, linetype = 2, color = "grey") + 
  labs(y = "-Log(p-value)", x = "Feature", title = "Significance level of attributes") + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank(), axis.text.x=element_blank()) + 
  scale_color_manual(values = cbPalette[2:3]) +
  geom_hline(yintercept = -log(0.05, 10), linetype = 4, color = "grey") 
t4
```

![](STIRexample_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


