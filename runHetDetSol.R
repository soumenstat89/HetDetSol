#####################################################################################################
#####   SINGLE-SEASON SCR ACCOUNTING FOR AUTOCORRELATED DETECTION: SIMULATION AND MCMC SCRIPT   #####
#####                           CONTINUOUS or CATEGORICAL EFFECT ON P0                          #####
#####################################################################################################

#' @title Function for simulation of SCR data set from SARE model amd to run MCMC
#'
#' @description
#' \code{runHetDetSol} returns the MCMC samples of monitored variables (myNimbleOutput), summary statistics of the parameters of ineterest (summary.output) and the simulation parameters (SIMULATION.PARAMETERS).
#' 
#' @param modelToFit A \code{Character} variable denoting the model to be fitted. The choices are: "SCR" (a basic single-season SCR model),
#'  "RE" (SCR model with independent random effects), "SARE" (SCR model with spatially autocorrelated random effects),
#'  "FM" (SCR model with two-group finite mixture), "FE" (SCR model with known true fixed effects).
#' @param ToAggregate A \code{LOGICAL} variable to indicate whether to aggregate the random effects.
#' @param sim.type \code{Character} variable denoting the type of simulation: "CON": to indicate continuous variation,
#'  "CAT" to indicate categorical variation in baseline detection probability.
#' @param WD A \code{Character} giving the working directory where the output will be saved.
#' @param N A \code{Numeric} variable denoting the population size.
#' @param M A \code{Numeric} variable denoting an upper bound of the population size (i.e., the size of the augmented SCR data set).
#' @param eta A \code{Numeric}  variable denoting the average baseline detection probability in SARE model to simulate SCR data.
#' @param sigma A \code{Numeric}  variable denoting the scale parameter of the half-normal detection function in SARE model to simulate SCR data.
#' @param phi A \code{Numeric}  variable controlling the level of autocorrelation in detection probability in SARE model to simulate SCR data.
#' @param niter A \code{Numeric} variable denoting number of MCMC iterations.
#' @param nburnin A \code{Numeric} variable denoting the burn-in period in MCMC computation.
#' @param nchains A \code{Numeric} variable denoting the number of chains in MCMC computation.
#' @param thin A \code{Numeric} variable denoting the thinning parameter in MCMC computation.
#' @param plot.check A \code{LOGICAL} variable to indicate whether the plots need to be saved in the working directory (WD) or not.

runHetDetSol <- function(
  modelToFit = "FE"
  , ToAggregate = F
  , sim.type = "CON" #"CAT"
  , WD = NA
  , N = 300
  , M = 500
  , eta = 0.3
  , sigma = 1.5
  , phi = 0.05
  , niter = 100
  , nburnin = 50
  , nchains = 2
  , thin = 1 
  , plot.check = FALSE){ ##--DO ALL
  
  require(rgeos)              # Import the geospatial analysis libraries
  require(rgdal)              # Import the spatial input/ouput libraries
  require(raster)             # Import the raster spatial package
  require(coda)               # Import the MCMC diagnostic tools
  require(nimble)             # Import the NIMBLE subroutines
  require(R.utils)            # Import the utility functions (some of the API of this package is experimental)
  require(abind)              # Import the library for manipulating multidimensional arrays
  require(nimbleSCR)
  
  
  if(is.na(WD)){WD <- getwd()}
  ts = format(Sys.time(), "%y%m%d_%H%M%S")
  if(!ToAggregate){ modelName = paste("HetDetSol_",modelToFit, "_",sim.type, "_", ts, sep = "")
  }else if(ToAggregate & !modelToFit %in% c("SCR", "FE")){modelName = paste("HetDetSol_",modelToFit, "_Aggregated_",sim.type,"_",  ts, sep = "")
  }else if(ToAggregate & modelToFit %in% c("SCR", "FE")){modelName = paste("HetDetSol_",modelToFit, "_",sim.type, "_", ts, sep = "")
  }
  if(!dir.exists(file.path(WD,modelName))){dir.create(file.path(WD, modelName))}
  
  
  ### SIMULATION PARAMETERS
  extent = 32    # habitat extent
  habRes = 2     # habitat resolution   
  buffer = 4.5   # buffer around the habitat 
  detRes = 1     # detector resolution 
  n.trials = 1   # number of samplign occasions
  
  if(plot.check){
    graphics.off()
    path <- file.path(WD,modelName, paste0("PLOTS", ".pdf"))
    pdf(file=path, width = 10, height = 7)
  }
  

  
  ## ----------------------------------------------------------------------------------------------
  ## ------ 1. SET-UP HABITAT AND DETECTORS -----
  ## ----------------------------------------------------------------------------------------------

  ## ------     1.1. GENERATE HABITAT -----
  coords <- matrix(c(buffer            , buffer ,
                     extent + buffer, buffer ,
                     extent + buffer, extent + buffer,
                     buffer            , extent + buffer,
                     buffer            , buffer ), ncol = 2, byrow = TRUE)
  
  P1 <- Polygon(coords)
  myStudyArea <- SpatialPolygons(list(Polygons(list(P1), ID = "a")),
                                 proj4string = CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  myHabitat <- MakeHabitat(poly = myStudyArea,
                           resolution = habRes,
                           buffer = buffer) 
  numHab <- sum(myHabitat$habitat.r[ ])
  habQualityDims <- dim(myHabitat$habitat.r)[2:1]
  
  # coordinates(myHabitat$habitat.r)
  
  ## ------     1.2. GENERATE DETECTORS -----
  myDetectors <- MakeSearchGrid( data = myStudyArea,
                                 resolution = detRes, 
                                 center = TRUE)                                        
  numDet <- nrow(myDetectors$detector.sp)
  
  if(plot.check){
    #--- PLOT DETECTORS
    plot(raster::union(myStudyArea, rgeos::gBuffer(myStudyArea, width = buffer)), col = "gray80")						
    plot(myStudyArea, col=grey(0.3), add=TRUE)
    plot(myDetectors$grid.poly, add=TRUE, lwd=1)
    plot(myDetectors$detector.sp, pch=19, cex=0.4, add=TRUE,
         col="cyan")
    axis(1)
    axis(2)
    
  }
  
  ## ------     1.3. SCALE DETECTORS & HABITAT UPPER/LOWER COORDINATES -----
  scaledDet <- UTMToGrid( data.sp = myDetectors$detector.sp,
                          grid.sp = myHabitat$habitat.sp)

  
  scaledHab <- UTMToGrid( data.sp = myHabitat$habitat.sp[myHabitat$habitat.r[]==1, ],
                          grid.sp = myHabitat$habitat.sp)$data.scaled.xy
  scaledUpperCoords <- scaledHab + 0.5 
  scaledLowerCoords <- scaledHab - 0.5 
  
  
  # CREATE A HABITAT GRID 
  habIDCells.mx <- myHabitat$IDCells.mx 
  habIDCells.mx[] <- 0
  scaledHabGridCenters <- scaledDet$grid.scaled.xy[myHabitat$habitat.r[]==1, ]
  for(ii in 1:nrow(scaledHabGridCenters)){
    habIDCells.mx[trunc(scaledHabGridCenters[ii,2])+1,
                  trunc(scaledHabGridCenters[ii,1])+1] <- ii
  }
  
  
  
  ## ------     1.4. CREATE CACHED OBJECTS FOR LESS RESTRICTION -----
  DetectorIndexLESS <- getLocalObjects(
    habitatMask = myHabitat$habitat.mx,
    coords = scaledDet$data.scaled.xy,
    dmax = 10/res(myHabitat$habitat.r)[1],
    resizeFactor = 1,
    plot.check = F
  )
  
  ## ----------------------------------------------------------------------------------------------
  ## ------ 2. SIMULATING SCR DATA SET -----
  ## ----------------------------------------------------------------------------------------------
  
  ## ------     2.1. SIMULATE AUTOCORRELATED COVARIATE -----
  distance.det <- gDistance( myDetectors$detector.sp,
                             myDetectors$detector.sp,
                             byid = T)
  
  V.mat <- exp(-phi * distance.det)

  logit.p0Vec <- rmnorm_chol( n = 1,
                              mean = rep(logit(eta), numDet),
                              cholesky = chol(V.mat),
                              prec_param = FALSE)
  
  if(sim.type == "CAT"){
    temp <- logit.p0Vec 
    ## Identify the cutoff value between the two categories so that "det.prop" % are
    ## in category one and 1-det.prop % are in category 2
    cutoff <- quantile(temp, 0.5)
    temp[logit.p0Vec <= cutoff] <- logit(0) # 1
    temp[logit.p0Vec > cutoff]  <- logit(eta) #2
    logit.p0Vec <- temp
    
  }
  
  p0Vec <- ilogit(logit.p0Vec)
  
  r.covs <- rasterFromXYZ(cbind(coordinates(myDetectors$detector.sp), logit.p0Vec))
  # r.covs <- rasterFromXYZ(cbind(coordinates(myDetectors$detector.sp), logit.p0Vec))
  w.movingwindow <- matrix(c(1,1,1,
                             1,0,1,
                             1,1,1), nrow = 3)
  MoransI <- Moran(r.covs, w.movingwindow)
  names(r.covs) <-  paste("Moran's Index",
                          round(MoransI, digits = 2), sep = "") 
  
  MyDensityRaster <- list(raster = r.covs,
                          distance = distance.det,
                          w.movingwindow = w.movingwindow,
                          # gamma = gamma,
                          phi = phi,
                          MoransI = MoransI)
  
  
  ##-- extract detector covs for a continuous cov
  myDetectors$detector.sp$detCov <- data.frame(raster::extract(MyDensityRaster$raster, myDetectors$detector.sp))[1] ##--[EM]
  
  ##--change name to detCov to be consistent over different runs
  colnames(myDetectors$detector.sp$detCov) <- "detCov"
  summary(myDetectors$detector.sp$detCov) ##--check
  dim(myDetectors$detector.sp$detCov)
  
  
  ## ------     2.2. SIMULATE INDIVIDUAL AC LOCATIONS -----
  ## Sample random (uniform) activity center locations
  ## sample ACs from the same habitat extent than the one used in the model 
  habIntensity <- rep(1, nrow(scaledLowerCoords))
  sxy <- matrix(NA, nrow = N, ncol = 2)
  for(ii in 1:N){
    sxy[ii, ]<- rbernppAC(n = 1,
                          lowerCoords = scaledLowerCoords,
                          upperCoords = scaledUpperCoords,
                          logIntensities = log(habIntensity),
                          logSumIntensity = log(sum(habIntensity)),
                          habitatGrid = habIDCells.mx,
                          numGridRows = nrow(myHabitat$IDCells.mx),
                          numGridCols = ncol(myHabitat$IDCells.mx))
    
  }#i
  
  # CONVERT TO SPATIAL POINTS
  simulated.ACS <- SpatialPoints(sxy)
  simulated.ACS$id <- 1:length(simulated.ACS)
  
  
  
  ## ------     2.3. SIMULATE DETECTION ARRAY  -----
  scaledDet.sp <- SpatialPoints(scaledDet$data.scaled.xy)
  scaledDet.sp$main.cell.id <- c(1:length(scaledDet.sp)) 
  scaledDet.sp$main.cell.x <- coordinates(scaledDet.sp)[ ,1] 
  scaledDet.sp$main.cell.y <- coordinates(scaledDet.sp)[ ,2] 
  
  scaledSigma <- sigma/habRes
  scaledDet.xy <- scaledDet$data.scaled.xy
  n.detectors <- nrow(myDetectors$detector.sp)
  Yarray.all <- array(0, c(N, n.detectors, n.trials))
  for(k in 1:n.trials){
    
    for(ii in 1:N){  
      #--- Calculate squared distance to all detectors
      d2 <- (scaledDet.xy[ ,1]-sxy[ii,1])^2 + (scaledDet.xy[ ,2]-sxy[ii,2])^2
      #--- Calculate detector-specific detection probability (half-normal detection function)
      p <- p0Vec * exp(-d2/(2*scaledSigma^2))
      #--- Sample binary individual detections 
      Yarray.all[ii,,k ] <- rbinom(n = length(p), size = 1, p)
    }#i 
  }#k
  
  y.all <- apply(Yarray.all, c(1,2), sum)
  detected <- apply(y.all,1, function(x) sum(x)>0)
  # cat("n.detected = ", sum(detected), "\n", sep = "")
  # cat("n.detections = ", sum(y.all), "\n", sep = "")
  # 
  detector.xyUnscaled <- coordinates(myDetectors$detector.sp)
  y <- y.all[detected,]
  
  ndetections.per.trap <- colSums(y)
  
  ## ------  2.4. To create some utility variables for aggregating random effects  -----
  if(ToAggregate & modelToFit %in% c("RE", "SARE", "FM")){  
    ntrapsPerZone <- 16
    B.list <- list()
    detectors.xy.zone.list <- list()
    distance.zone.list <- list()
    c3 <- 0
    
    nZones <- n.detectors / ntrapsPerZone
    trapsInZones <- fun.trapsInZones(n.detectors=n.detectors, ntrapsPerZone=ntrapsPerZone)
    
    detectors.xy.zone <- matrix(NA, nZones, 2)
    dimnames(detectors.xy.zone) <-  dimnames(detector.xyUnscaled)
    ndetections.per.zone <- rep(NA, nZones)
    for(j in 1:nZones){
      detectors.xy.zone[j,] <- apply(detector.xyUnscaled[trapsInZones[,j],],2,mean)
      ndetections.per.zone[j] <- sum(ndetections.per.trap[trapsInZones[,j]])
      
    }
    
    # Tranformation matrix from Zonal-Random-Effect to Detector-Specific-Random-Effect
    B <- matrix(0, n.detectors, nZones)
    for(j in 1:nZones){
      B[trapsInZones[,j] ,j] <- 1
    }
    
    c3<- c3+1
    B.list[[c3]] <- B
    detectors.xy.zone.list[[c3]] <- detectors.xy.zone
    
    distance.zone <- sqrt(e2dist(detectors.xy.zone,detectors.xy.zone))
    distance.zone.list[[c3]] <- distance.zone
    
    
    names(B.list) <- c("B.16TrapsPerZone")
    names(detectors.xy.zone.list) <- c("detectors.xy.zone.16TrapsPerZone")
    names(distance.zone.list) <- c("distance.zone.16TrapsPerZone")
    
  }#ToAggregate
  
  ## ------     2.5. AUGMENT DATA -----
  
  ## augmentation factor that ensures that the total number of individuals 
  ## (alive + available) is always the same, regardless of the simulation
  this.AF <- M/sum(detected)-1
  
  sum(detected) * (1+ this.AF) ==  M ##--check
  
  y <- MakeAugmentation(y = y, aug.factor = this.AF, replace.value = 0)
  z <- MakeAugmentation( y = rep(1, sum(detected)),
                         aug.factor = this.AF, replace.value = NA)
  length(z)#dim(z)
  
  
  ## ------     2.6. TRANSFORM Y TO SPARSE MATRICES -----
  
  SparseY <- getSparseY(y, nMaxTraps = 100)
  
  ## ------     2.7. SET THE INPUT FOR NIMBLE -----
  
  
  cat("Model: ", modelToFit, "\n", sep="")
  if(ToAggregate & modelToFit %in% c("RE", "SARE", "FM")){ 
    cat("Aggregating random effects...", "\n", sep="")
  }
  nimConstants <- list(
    n.detectors = dim(y)[2], 
    n.individuals = dim(y)[1],
    numHabWindows = numHab,
    x.max = habQualityDims[1],
    y.max =  habQualityDims[2],
    y.maxDet = dim(DetectorIndexLESS$habitatGrid)[1], 
    x.maxDet = dim(DetectorIndexLESS$habitatGrid)[2], 
    ResizeFactor = DetectorIndexLESS$resizeFactor, 
    dimYbinded = SparseY$lengthYCombined, 
    maxNBDets = DetectorIndexLESS$numLocalIndicesMax, 
    n.cells = dim(DetectorIndexLESS$localIndices)[1])
  
  
  nimData <- list(
    
    z = z,
    detector.xy = scaledDet$data.scaled.xy,
    yBinded = SparseY$yCombined[,,1],
    trials = rep(n.trials, nrow(scaledDet$data.scaled.xy)),  
    localTrapsIndices = DetectorIndexLESS$localIndices, 
    numLocalTraps = DetectorIndexLESS$numLocalIndices, 
    habitatGridDet = DetectorIndexLESS$habitatGrid,
    habitatGrid = habIDCells.mx, 
    lowerHabCoords = scaledLowerCoords,
    upperHabCoords = scaledUpperCoords,
    logHabIntensity = log(habIntensity),
    logSumHabIntensity = log(sum(habIntensity)),
    distance = distance.det)
  
  
  sxy.init <- MakeInitXY( y = y,
                          habitat.mx = myHabitat$habitat.mx,
                          detector.xy = scaledDet$data.scaled.xy,
                          IDCells.mx = myHabitat$IDCells.mx,
                          grid.xy = scaledDet$grid.scaled.xy)
  
  ## Initialise z values
  z.init <- ifelse(!is.na(z), NA, 1)
  z.init[!is.na(z.init)] <- rbinom(sum(!is.na(z.init)), size = 1, prob = 0.5)
  
  
  
  log.phi.init <-  log(runif(1, 0.1, 2))
  logit.eta.init <- logit(runif(1, 0.1, 0.6))
  V.mat <- exp(-exp(log.phi.init) * distance.det)
  W.init <- rmnorm_chol(n=1,
                        mean=rep(logit.eta.init, n.detectors),
                        cholesky = chol(V.mat),
                        prec_param = FALSE) 
  nimInits.orig <- list(
    z = z.init,    
    psi = runif(1,0,1),
    sxy = sxy.init, #[ , ,1],
    logit.eta = logit.eta.init, 
    sigma = runif(1,0,10), 
    log.phi =log.phi.init,
    W = W.init)
  
  if(modelToFit=="SCR"){ ##-- SCR model 
    modelCode <- modelCode_SCR
    
    remove.index <- c("distance")
    indices <- which(names(nimData)%in%remove.index)
    if(length(indices)>0) nimData <- nimData[-indices]
    
    nimInits <- nimInits.orig
    nimInits$logit.p0 <- logit(runif(1, 0.1, 0.6))
    remove.index <- c("W", "log.phi", "logit.eta")
    indices <- which(names(nimInits)%in%remove.index)
    if(length(indices)>0) nimInits <- nimInits[-indices]
    
    nimParams <- c("p0", "logit.p0", "sigma","p0Vec",
                   "psi","N") # , "z","sxy"
    
  }else if(!ToAggregate & modelToFit == "SARE"){ ##--SARE
    modelCode <- modelCode_SARE
    
    nimInits <- nimInits.orig
    
    nimParams <- c("log.phi","phi", "logit.eta", "eta", "sigma","p0Vec", # "W",
                   "psi", "z","sxy","N") # , "z","sxy"
  }else if(!ToAggregate & modelToFit=="RE"){ ##-- RE model 
    modelCode <- modelCode_RE
    
    remove.index <- c("distance") 
    indices <- which(names(nimData)%in%remove.index)
    if(length(indices) >0) nimData <- nimData[-indices]  
    
    nimInits <- nimInits.orig
    nimInits$sigma.w <- rgamma(1, shape = 1, rate = 1)
    remove.index <- c("log.phi") 
    indices <- which(names(nimInits)%in%remove.index)
    if(length(indices) >0) nimInits <- nimInits[-indices]  
    
    nimParams <- c("logit.eta", "eta", "sigma", "p0Vec", #"W",
                   "sigma.w",  "psi","N") # , "z","sxy"
    
    
  }else if(!ToAggregate & modelToFit=="FM"){ ##-- FM model 
    modelCode <- modelCode_FM
    
    n.groups <- 2    ## fixed number of groups
    nimConstants$n.groups <- n.groups
    
    
    nimData$one <- rep(1,n.groups-1)
    remove.index <- c("distance") 
    indices <- which(names(nimData)%in%remove.index)
    if(length(indices) >0) nimData <- nimData[-indices]  
    
    y.all <- apply(Yarray.all, c(1,2), sum)
    ndetections.per.trap <- colSums(y.all)
    trapsWithDetections <- ndetections.per.trap > 0
    u.init <- rep(0, n.detectors)
    u.init[trapsWithDetections] <- 1
    nimInits <- nimInits.orig
    nimInits$u <- u.init
    nimInits$pi <- sum(u.init == 1)/n.detectors
    nimInits$eta <- c(runif(1, 0, 0.2), runif(1, 0.21, 0.5)) 
    remove.index <- c("logit.eta", "log.phi", "W") 
    indices <- which(names(nimInits)%in%remove.index)
    if(length(indices) >0) nimInits <- nimInits[-indices]  
    
    nimParams <- c( "eta", "sigma", "psi","N", "pi", "u", "p0Vec") # , "z","sxy"
    
    
  }else if(ToAggregate & modelToFit == "SARE"){ ##--SARE with aggregated random effects
    modelCode <-  modelCode_SARE.agg
    
    nimConstants$nZones <- nZones
    
    nimData$B <- B.list[[paste0("B.",ntrapsPerZone, "TrapsPerZone")]]
    nimData$distance.zone <- distance.zone.list[[paste0("distance.zone.",ntrapsPerZone, "TrapsPerZone")]]
    remove.index <- c("distance")
    indices <- which(names(nimData)%in%remove.index)
    if(length(indices) >0) nimData <- nimData[-indices]  
    
    nimInits <- nimInits.orig
    V.mat <- exp(-exp(nimInits$log.phi) * nimData$distance.zone)
    nimInits$Wzone <- rmnorm_chol(n=1,
                                  mean=rep(nimInits$logit.eta, nZones),
                                  cholesky = chol(V.mat),
                                  prec_param = FALSE) 
    
    remove.index <- c("W") 
    indices <- which(names(nimInits)%in%remove.index)
    if(length(indices) >0) nimInits <- nimInits[-indices] 
    
    nimParams <- c("log.phi","phi", "logit.eta", "eta", "sigma", "p0Vec", # "Wzone",
                   "psi","N")# , "z","sxy"
    
  }else if(ToAggregate & modelToFit=="RE"){ ##--RE with aggregated random effects
    modelCode <- modelCode_RE.agg
    
    nimConstants$nZones <- nZones
    
    nimData$B <- B.list[[paste0("B.",ntrapsPerZone, "TrapsPerZone")]]
    distance.zone <- distance.zone.list[[paste0("distance.zone.",ntrapsPerZone, "TrapsPerZone")]]
    remove.index <- c("distance")
    indices <- which(names(nimData)%in%remove.index)
    if(length(indices) >0) nimData <- nimData[-indices]  
    
    
    nimInits <- nimInits.orig
    V.mat <- exp(-1 * distance.zone)
    nimInits$Wzone <- rmnorm_chol(n=1,
                                  mean=rep(nimInits$logit.eta, nZones),
                                  cholesky = chol(V.mat),
                                  prec_param = FALSE) 
    nimInits$sigma.w <- rgamma(1, shape = 1, rate = 1)
    remove.index <- c("W","log.phi") 
    indices <- which(names(nimInits)%in%remove.index)
    if(length(indices) >0) nimInits <- nimInits[-indices]
    
    nimParams <- c("logit.eta", "eta", "sigma", "p0Vec", #"Wzone",
                   "sigma.w",  "psi", "N")# , "z","sxy"
    
    
  }else if(ToAggregate & modelToFit=="FM"){ ##--FM with aggregated random effects
    modelCode <-  modelCode_FM.agg
    
    n.groups <- 2    ## fixed number of groups
    nimConstants$nZones <- nZones
    nimConstants$n.groups <- n.groups
    
    nimData$B <- B.list[[paste0("B.",ntrapsPerZone, "TrapsPerZone")]]
    nimData$one <- rep(1,n.groups-1)
    remove.index <- c("distance") 
    indices <- which(names(nimData)%in%remove.index)
    if(length(indices) >0) nimData <- nimData[-indices]  
    
    zonesWithDetections <- ndetections.per.zone > 0
    nimInits <- nimInits.orig
    u.init <- rep(0, nZones)
    u.init[zonesWithDetections] <- 1
    nimInits$u <- u.init
    nimInits$pi <- sum(u.init == 1)/nZones
    nimInits$eta <- c(runif(1, 0, 0.2), runif(1, 0.21, 0.5)) 
    
    remove.index <- c("logit.eta", "log.phi", "W") 
    indices <- which(names(nimInits)%in%remove.index)
    if(length(indices) >0) nimInits <- nimInits[-indices]  
    
    nimParams <- c("eta", "sigma","psi", "N", "pi", "u", "p0Vec")#, "p0Zone")# , "z","sxy"
    
    
  }else if(modelToFit == "FE"){ ##-- FE model 
    
    modelCode <- modelCode_FE
    
    nimData$W <- c(myDetectors$detector.sp$detCov[,1])
    
    nimInits <- nimInits.orig
    
    remove.index <- c("W")
    indices <- which(names(nimInits)%in%remove.index)
    if(length(indices)>0) nimInits <- nimInits[-indices]
    
    nimParams <- c("log.phi","phi", "logit.eta", "eta", "sigma", "p0Vec",
                   "psi", "N")# , "z","sxy"
    
  }##
  
  
  
  
  ## -----------------------------------------------------------------------------
  ## ------ 3. FIT NIMBLE MODEL -----
  ## -----------------------------------------------------------------------------
  
  ## ------     3.1. BUILDING NIMBLE MODEL -----
  
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        data = nimData,
                        inits = nimInits,
                        check = FALSE,
                        calculate = FALSE)
  cat("Joint density at initial value:\n", sep = "")
  print(model$calculate())
  model$initializeInfo()
  
  ## ------     3.2. COMPILING NIMBLE MODEL -----
  cmodel <- compileNimble(model)
  print(cmodel$calculate())
  
  ## ------     3.3. ASSIGN MCMC SAMPLERS TO THE NODES -----
  MCMCconf <- configureMCMC( model = model,
                             monitors = nimParams,
                             control = list(reflective = TRUE),
                             thin = 1) 
  
  if(!ToAggregate & modelToFit %in% c("FM")){ #, "FM.agg")){
    MCMCconf$removeSamplers("u")
    
    MCMCconf$addSampler(target = "u", type = "binaryFEWatAtime", numIndicesToChange = 4) ### CORRECT as I want to add binaryFEWatAtime sampler to vector u[1:n.detector]
    # # # retrieve the current ordering of sampler execution
    ordering <- MCMCconf$getSamplerExecutionOrder()
    # ordering
    len <- length(ordering)
    MCMCconf$setSamplerExecutionOrder(c(rep(ordering[1], 10), # eta[1]
                                        rep(ordering[2], 10), # eta[2]
                                        rep(ordering[3], 10), # pi
                                        ordering[4:(len-1)],
                                        rep(ordering[len], 10)))#, # u[1:n.detectors]
  }#if(modelToFit %in% c("FM")){
  
  
  ## ------     3.4. BUILDING MCMC CONFIGURATION -----
  MCMC <- buildMCMC(MCMCconf)
  
  ## ------     3.5. COMPILE MCMC -----
  
  cMCMC <- compileNimble( MCMC,
                          project = model,
                          resetFunctions = TRUE)
  
  ## ------     3.6. RUN THE MCMC -----
  
  
  Runtime <- system.time(myNimbleOutput <- runMCMC( cMCMC,
                                                    niter = niter,
                                                    thin = thin, 
                                                    nburnin = nburnin,
                                                    nchains = nchains,  
                                                    samplesAsCodaMCMC = TRUE))
  
 
  cat("Time taken: ", Runtime[3]/60, " minutes", "\n\n", sep = "")
  
  
  ## -----------------------------------------------------------------------------
  ## ------ 4. PROCESS MCMC OUTPUT -----
  ## -----------------------------------------------------------------------------
  
  ## ------     4.1. PLOT THE OUTPUT  -----
  paramnames <- c('N', 'sigma')
  if(is.list(myNimbleOutput)){
    myNimbleOutput2 <-  as.mcmc.list(lapply(myNimbleOutput, function(x){
      out <-  mcmc(x[,paramnames])
      out[,"sigma"] <- out[,"sigma"]*habRes  ##-- Unscaling the sigma
      return(out)
    }))
  }
  if(!is.list(myNimbleOutput)){
    myNimbleOutput2 <-  mcmc(myNimbleOutput[,paramnames])
    myNimbleOutput2[,"sigma"] <- myNimbleOutput2[,"sigma"]*habRes  ##-- Unscaling the sigma
    
  }
  
  if(plot.check){
    #--- MCMC TRACEPLOTS AND DENSITY
    # basicMCMCplots::chainsPlot(myNimbleOutput2, line = c(N, sigma))
    
    sim.values <- c(N, sigma)
    par(mfrow=c(2,2))
    ii <- 1
    for(ii in 1:length(paramnames)){
      cols <- c("darkturquoise", "deeppink", "darkorange", "blueviolet",  "darkcyan","red") #, "royalblue4")
      traceplot(myNimbleOutput2[ ,paramnames[ii]], col = cols)
      abline(h = sim.values[ii], col = "black", lwd = 3, lty = 2)
      plot(density(unlist(myNimbleOutput2[ ,paramnames[ii]])), main = paramnames[ii])
      abline(v = sim.values[ii], col = "red", lwd = 2)
    }
    
  }#plot check 
  
  ## ------     4.2. SUMMARY OF THE OUTPUT  -----
  # myNimbleOutput <- as.mcmc.list(myNimbleOutput) 
  if(is.list(myNimbleOutput2)){posteriorSamples <- do.call(rbind, myNimbleOutput2)} 
  if(!is.list(myNimbleOutput2)){posteriorSamples <- as.matrix(myNimbleOutput2)} 
  
  summary.output <- do.call(rbind, lapply(paramnames, function(x){
    this.posterior <- posteriorSamples[, x]
    out <- list(MEAN = mean(this.posterior),
                CV = sd(this.posterior)/mean(this.posterior),
                LCI = quantile(this.posterior, 0.025),
                UCI = quantile(this.posterior, 0.975))
    return(unlist(out))
  }))
  dimnames(summary.output)[[1]] <- paramnames
  
  # cat("No. of detections = ", sum(y.all), "\n", sep = "")
  
  cat('No. of detected individuals = ',sum(detected), '\n\n', sep = '')
  
  cat('Posterior summary output:', '\n\n', sep = '')
  print(summary.output)
  
  
  
  ### SIMULATION PARAMETERS
  SIMULATION.PARAMETERS <- list(
    # habitat extent
    HABITAT.extent = extent,
    # habitat resolution
    HABITAT.resolution = habRes,
    # buffer around the habitat
    HABITAT.buffer = buffer,
    # detector resolution
    DETECTOR.resolution = detRes,
    # population size
    N = N,
    # size of the augmented population
    M = M,
    # baseline detection probability
    eta = eta,
    # the scale parameter
    sigma = sigma,
    # autocorrelation parameter
    phi =  phi,
    # number of occasion
    n.trials = n.trials,
    # continuous or categorical variation
    sim.type = sim.type,
    # Whether to aggregate random effects
    ToAggregate = ToAggregate
  )
  
  
  ## ------     4.3. PLOT THE TRUE ANS ESTIMATED p0 SURFACE -----
  
  # myNimbleOutput <- as.mcmc.list(myNimbleOutput) 
  if(is.list(myNimbleOutput)){
    p0Vec.chain <- do.call(rbind, lapply(myNimbleOutput,function(x){
      p0Vec.names <- grep("p0Vec", colnames(x), value = T)
      return(x[, p0Vec.names])
    }))
  }else if(!is.list(myNimbleOutput)){
    p0Vec.chain <- as.matrix(myNimbleOutput)
    p0Vec.chain <- p0Vec.chain[, grep("p0Vec", colnames(p0Vec.chain), value = T)] 
  } 
  p0Vec.est <- colMeans(p0Vec.chain)
  
  if(plot.check){
    ##--- PLOT THE BASELINE DETECTION PROBABILITY SURFACE
    n.col <- 50 
    r.col1 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(5, "OrRd"))(n.col)
    col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(5, "GnBu"))(2)
    colfunc <- colorRampPalette(c(col[2],r.col1[1])) 
    r.col2 <- colfunc(100-n.col)
    r.col <- c(r.col2, r.col1)
    
    par(mfrow=c(1,2), mar=c(2,2,2,2), oma=c(2,2,1,3))
    
    range.p <- range(c(p0Vec, p0Vec.est))
    if(is.infinite(range.p[1])) range.p[1] <- -10 # -10^5
    cuts.p <- seq(range.p[1], range.p[2], length.out = 100)
    
    p0Vec.r <- rasterFromXYZ(cbind(coordinates(myDetectors$detector.sp), p0Vec))
    plot(p0Vec.r,
         axes = FALSE,
         box = FALSE,
         breaks = cuts.p,
         col = r.col,
         bty = "n",
         main = paste0("TRUE p0 Surface (Moran's.I = ", round(Moran(p0Vec.r), 2),")"),
         legend = F  
    ) #
    
    p0Vec.est.r <- rasterFromXYZ(cbind(coordinates(myDetectors$detector.sp), p0Vec.est))
    plot(p0Vec.est.r,
         axes = FALSE,
         box = FALSE,
         breaks = cuts.p,
         col = r.col,
         bty = "n",
         main = paste0("MEAN p0 Surface (Moran's.I = ", round(Moran(p0Vec.est.r), 2),")"),
         legend = F  #T,
         
    ) #
    
    # # plot.new()
    # par(xpd = TRUE)
    
    plot(p0Vec.r, breaks = cuts.p, col = r.col,
         legend.width = 2, legend.only = T,
         axis.args = list(at =     round(seq(range.p[1], range.p[2], length.out = 6), digits = 1),
                          labels = round(seq(range.p[1], range.p[2], length.out = 6), digits = 1),
                          cex.axis = 1),
         legend.args = list(text = '', # 'p0Vec',
                            line = -1, side = 4, font = 2, cex = 1))
    
  
    
  }#plot check
  if(plot.check){
    graphics.off()
  }
  
  ## ------     4.3. SAVE THE OUTPUT -----
  
  file.name <- paste( "NimbleOutFile_fit", modelToFit, ".RData", sep = "")
  this.path <- file.path(WD, modelName, file.name)
  
  if(ToAggregate & modelToFit %in% c("RE", "SARE", "FM")){
    save(myNimbleOutput, summary.output, 
         Runtime, niter, nburnin, nchains,
         modelToFit,paramnames,
         nimData, nimConstants, nimInits, nimParams,
         modelCode, MyDensityRaster, myDetectors, SIMULATION.PARAMETERS,
         simulated.ACS, Yarray.all, detector.xyUnscaled,
         B.list, detectors.xy.zone.list, distance.zone.list,
         file=this.path)
  }else{save(myNimbleOutput, summary.output, 
             Runtime, niter, nburnin, nchains,
             modelToFit,paramnames,
             nimData, nimConstants, nimInits, nimParams,
             modelCode, MyDensityRaster, myDetectors, SIMULATION.PARAMETERS,
             simulated.ACS, Yarray.all, detector.xyUnscaled,
             # B.list, detectors.xy.zone.list, distance.zone.list,
             file=this.path)
  }
  
}#runHetDetSol

