#=========================================================================#
# 
# -------------------- SIMULATE N.REP DATASETS ---------------------------#
# 
# -------------------[with attribute interactions]------------------------#
# 
#=========================================================================#

for (repi in (1:n.rep)){
  
  if (letsSimulate == TRUE){
    if (is.main){ # make this main effect sim if TRUE
      sim.data <- createSimulation(num.samples = n.samp,
                                   num.variables = num.attr,
                                   pct.signals = pct.signals,
                                   pct.train = 1 / 2,
                                   pct.holdout = 1 / 2,
                                   # pct.validation = 1 /3,
                                   bias = bias,
                                   sim.type = "mainEffect",
                                   verbose = FALSE)
    } else { # interaction simulation
      sim.data <- createSimulation(num.samples = n.samp,
                                   num.variables = num.attr,
                                   pct.signals = pct.signals,
                                   pct.train = 1 / 2,
                                   pct.holdout = 1 / 2,
                                   bias = bias,
                                   # pct.validation = 1 /3,
                                   #meanExpression = 7, A = NULL, randSdNoise = 1, sdNoise = 0.4,
                                   sim.type = "interactionErdos", # interactionScalefree
                                   verbose = FALSE)
    }
    dat <- rbind(sim.data$train, sim.data$holdout)
  } else {
    dat <- read.csv("simulatedInte.csv")
  }
  
  dat[, class.lab] <- as.factor(dat[, class.lab]) 
  pheno.class <- dat[, class.lab]
  predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
  attr.names <- colnames(predictors.mat)
  num.samp <- nrow(dat)
  neighbor.idx <- find.neighbors(predictors.mat, method = RF.method, pheno.class, k = k)
  
  if (writeData == TRUE){
    write.csv(dat, file = paste("data/ARF_compare", repi, RF.method, bias, "bias", 
                                ifelse(RF.method == "relieff", k, "No"), "k", 
                                pct.signals, "pct.signals",
                                num.attr, "num.attr", num.samp, "num.samp.csv", sep = "_"))
  }
}

