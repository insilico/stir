
# STIR (STatistical Inference Relief)

Trang T. Le, Ryan J. Urbanowicz, Jason H. Moore, B. A McKinney. “STatistical Inference Relief (STIR) feature selection,” Bioinformatics. 18 September 2018. [https://doi.org/10.1093/bioinformatics/bty788](https://doi.org/10.1093/bioinformatics/bty788)

[http://insilico.utulsa.edu/software/STIR](http://insilico.utulsa.edu/software/STIR)

Example STIR usage and output: [STIRexample.md](https://github.com/insilico/STIR/blob/master/inst/example/STIRexample.md).

### To install:

    >library(devtools)
    
    >install_github("insilico/stir")  # you can use build_vignettes = TRUE but slows down install

    >library(stir)
    
    >data(package="stir")
    
    # >vignette("STIRsimulated") # if you build_vignettes
    
    # >vignette("STIRmdd")

    
 ### Examples

[Simulated Data Example with privateEC Simulation](https://github.com/insilico/STIR/blob/master/inst/example/STIRexample.md) 

[RNA-Seq Example](https://github.com/insilico/STIR/blob/master/vignettes/STIRmdd.Rmd) 

### Abstract

#### Motivation

Relief is a family of machine learning algorithms that uses nearest-neighbors to select features whose association with an outcome may be due to epistasis or statistical interactions with other features in high-dimensional data. Relief-based estimators are non-parametric in the statistical sense that they do not have a parameterized model with an underlying probability distribution for the estimator, making it difficult to determine the statistical significance of Relief-based attribute estimates. Thus, a statistical inferential formalism is needed to avoid imposing arbitrary thresholds to select the most important features. 

#### Methods

We reconceptualize the Relief-based feature selection algorithm to create a new family of STatistical Inference Relief (STIR) estimators that retains the ability to identify interactions while incorporating sample variance of the nearest neighbor distances into the attribute importance estimation. This variance permits the calculation of statistical significance of features and adjustment for multiple testing of Relief-based scores. Specifically, we develop a pseudo t-test version of Relief-based algorithms for case-control data.  

#### Results

We demonstrate the statistical power and control of type I error of the STIR family of feature selection methods on a panel of simulated data that exhibits properties reflected in real gene expression data, including main effects and network interaction effects. We compare the performance of STIR when the adaptive radius method is used as the nearest neighbor constructor with STIR when the fixed-$k$ nearest neighbor constructor is used. We apply STIR to real RNA-Seq data from a study of major depressive disorder and discuss STIR's straightforward extension to genome-wide association studies.

#### Websites
[http://insilico.utulsa.edu/software/STIR](http://insilico.utulsa.edu/software/STIR)
[https://github.com/insilico](https://github.com/insilico)

#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
[ttle@pennmedicine.upenn.edu](ttle@pennmedicine.upenn.edu)
