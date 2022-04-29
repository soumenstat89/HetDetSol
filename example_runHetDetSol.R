## -----------------------------------------------------------------------------
# {##--DO ALL
  
  rm(list=ls())
  gc()
  
  ## source working directories 
  # doc.path <- "C:/Users/admin/OneDrive - Norwegian University of Life Sciences/Documents" ##--[SD]
  code.folder<- "E:/RovQuantSD_Codes/HetDetSol/MS_FINAL"
  # getwd()
  setwd(code.folder)
  
  source("HetDetSol_functions.R")
  # source("runHetDetSol.R")
  
  # mymodelFit <- c("SCR", "RE", "SARE", "FM", "FE")
  
  model.set <- c("SCR", "RE", "SARE", "FM", "FE")
  
  # modelToFit = "SCR"
  # ToAggregate = F
  # sim.type = "CAT" #"CAT"
  # WD = code.folder
  # N = 300
  # M = 500
  # eta = 0.6
  # sigma = 1.5
  # phi = 0.05
  # niter = 200
  # nburnin = 100
  # nchains = 1
  # thin = 1
  # plot.check = T
  
  # source("HetDetSol_utility_functions.R")
  
  runHetDetSol(
    modelToFit = "FE" #  "SARE"
    , ToAggregate = F
    , sim.type = "CON" #"CON" #"CAT" #
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
  
 