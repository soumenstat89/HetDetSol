
  rm(list=ls())
  gc()
  
  ## source working directories 
  code.folder<-  getwd()
  setwd(code.folder)

  ## source the HetDetSol functions
  source("HetDetSol_functions.R")
  model.set <- c("SCR", "RE", "SARE", "FM", "FE")
    
  runHetDetSol(
    modelToFit = "FE"
    , ToAggregate = F
    , sim.type = "CON" 
    , WD = code.folder
    , N = 300
    , M = 500
    , eta = 0.6
    , sigma = 1.5
    , phi = 0.05
    , niter = 200
    , nburnin = 100
    , nchains = 1
    , thin = 1
    , plot.check = T)
  
 
