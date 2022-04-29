library(nimble)                    # Import the NIMBLE subroutines
library(rgeos)                     # Import the geospatial analysis libraries
library(rgdal)                     # Import the spatial input/ouput libraries
library(raster)                    # Import the raster spatial package
library(coda)                      # Import the MCMC diagnostic tools
library(nimble)                    # Import the NIMBLE subroutines
library(nimbleSCR)
library(abind)                     # Import the library for manipulating multidimensional arrays
library(basicMCMCplots)


# require(rgeos)              # Import the geospatial analysis libraries
# require(rgdal)              # Import the spatial input/ouput libraries
# require(raster)             # Import the raster spatial package
# require(coda)               # Import the MCMC diagnostic tools
# require(nimble)             # Import the NIMBLE subroutines
# # require(ggplot2)            # Import the graphical libraries
# require(R.utils)            # Import the utility functions (some of the API of this package is experimental)
# require(abind)              # Import the library for manipulating multidimensional arrays
# require(nimbleSCR)
# # require(spatstat)
# # require(maptools)
# # require(parallel)
# # require(spatstat.geom)
# library(basicMCMCplots)


###############################################################
##################### MakeHabitat ###############################
###############################################################
#' @title MakeHabitat: a function to create a list of objects to define the habitat
#'
#'
#'
#' @description
#' \code{MakeHabitat} returns a \code{List} of objects necessary to define the habitat in SCR (can be slow for large habitat...)
#'
#' @param poly A \code{SpatialPolygon} object with the study area
#'	@param resolution A \code{Numerical} with the grid cell size of the habitat raster.
#' @param buffer A \code{Numeric} with the size of the buffer area around the focal study area (in meters in poly is UTM).
#'
#' 
#' @return A \code{List} object with the coordinates and attributes of the habitat
#'  @return habitat.sp: A \code{SpatialPointsDataFrame} object with the coordinates of the center of cells of the habitat
#'  @return habitat.clip.sp: A \code{SpatialPointsDataFrame} object with the coordinates of the center of cells of the habitat cliped to the suitable habitat.
#'  @return habitat.xy: A \code{DataFrame} object with the x and y coordinates of the habitat cells.
#'  @return IDCells.mx: A \code{matrix} with the ID of each cell of the raster
#'  @return habitat.r: A \code{raster} with the suitable habitat (1) and non suitable (0)
#'  @return habitat.mx: A \code{matrix} with the suitable habitat (1) and non suitable (0)
#'  @return resolution: A \code{numeric} with the resolution used for the raster 
#'  @return buffered.habitat.poly: A \code{SpatialPolygons} that includes the buffer. 
#'  @return habitat.poly: A \code{SpatialPolygons} with the original polygon used. 
#'  @return habitat.index: A \code{vector} with ID cell of the raster that fall within the suitable habitat. 



MakeHabitat <- function(poly, 			               ## Polygon of the study area (sp object)
                        resolution, 	               ## Resolution of the grid cell (in meters)
                        buffer = NULL)#, 					## Add a buffer size (in meters)
  # polygon.clip = NULL, 		   ## Add a polygon that would clip over the buffer (for the sea for example )
  # plot.check = TRUE, 			   ## if TRUE some graphics will show up to check if it looks right )
  # fasterize = FALSE,            ## Increase spee, requires library(fasterize)    
  # CoverToKeepHabitat = NULL)
{ 
  
  # poly <- myStudyArea
  # resolution <- myVars$HABITAT$resolution
  
  ## ----- Store the original study area polygon -----
  polygon.orig <- poly
  
  ## ----- Create the buffer zone around the study area -----
  buff <- rgeos::gBuffer(poly, width = buffer)
  poly <- raster::union(poly, buff)
  
  ## ----- Create Habitat raster (study area + buffer)-----   
  r <- raster(extent(poly))                                         ## Rasterise polygon
  res(r) <- resolution  ## Change resolution
  
  polygon.r <- rasterize(as(poly, "SpatialLines"), r, field=1)      ## rasterise as lines and then polygons so all cells touched by the polygon are covered 
  polygon.r <- rasterize(poly, polygon.r, update=T, field=1)        ## Add field=1 so its 1 for study area and NA for the rest.
  
  ## ----- Create Habitat matrix (habitat : 1 and non-habitat: 0) -----
  habitat.r <- polygon.r
  habitat.r[is.na(habitat.r)] <- 0                                     ## Give 0 values to raster cells outside the study area
  habitat.mx <- as.matrix(habitat.r)                                   ## Convert to matrix
  
  # ## ----- Give unique IDs to cells ----- 
  IDCells.r <- polygon.r
  IDCells.r[] <- 1:length(IDCells.r)							         ## Cell ID starts from the top left corner
  IDCells.mx <- as.matrix(IDCells.r)							               ## Convert to matrix
  
  ## ----- Obtain xy coordinates of cells -----   
  habitat.xy <- xyFromCell(polygon.r, 1:ncell(polygon.r))
  dimnames(habitat.xy) <- list(1:length(habitat.xy[,"x"]), c("x","y"))
  habitat.sp <- SpatialPointsDataFrame(data.frame(habitat.xy[,c("x","y")]), data=data.frame(habitat.xy), proj4string=CRS(projection(poly)))
  habitat.index <- which(!is.na(over(as(habitat.sp, "SpatialPoints"), as(poly,"SpatialPolygons"))))
  habitat.clip.sp <- habitat.sp[!is.na(over(as(habitat.sp, "SpatialPoints"), as(poly, "SpatialPolygons"))), ]
  
  ## -- Obtain lower and upper cell coordinates
  lower.hab.sp <- data.frame(coordinates(habitat.sp) - resolution/2)
  upper.hab.sp <- data.frame(coordinates(habitat.sp) + resolution/2)
  coordinates(upper.hab.sp) <- upper.hab.sp
  coordinates(lower.hab.sp) <- lower.hab.sp
  proj4string(lower.hab.sp) <- CRS(projection(poly))
  proj4string(upper.hab.sp) <- CRS(projection(poly))
  
  myHabitat <- list(habitat.sp = habitat.sp, #
                    habitat.clip.sp = habitat.clip.sp,
                    habitat.xy = habitat.xy,
                    IDCells.mx = IDCells.mx,#
                    habitat.r = habitat.r, #
                    habitat.mx = habitat.mx,
                    resolution = resolution,
                    buffered.habitat.poly = poly,#
                    habitat.poly = polygon.orig,#
                    habitat.index = habitat.index,
                    upper.hab.sp = upper.hab.sp, ##
                    lower.hab.sp = lower.hab.sp ##
  )
  return(myHabitat)
  
}

###############################################################
##################### MakeSearchGrid ###############################
###############################################################

#' @title Function to create detectors
#' #'
#' @description
#' \code{MakeSearchGrid} creates  detectors within a defined area (data) and resolution between detectors. 
#' Division precises if subspatial division should be performed (Following PAB model).  
#' 
#' @param data A \code{sp} object. can be a \code{SpatialPolygons} or \code{SpatialPoints} or \code{SpatialLines} or \code{raster}. Spatial Lines objects takes time to compute
#' @param resolution Numeric variable denoting the size of grid cells in units of (\code{polygon}).
#' @param center Logical variable denoting whether nodes of the resulting search grid are to be centered within their respective subplots.


MakeSearchGrid <- function(data,
                           resolution,
                           # div = 1,
                           center = T#,
                           # plot = TRUE,
                           # fasterize = FALSE
){
  
  ### ====  GENERATE DETECTORS ====
  
  # data <- myStudyArea
  # resolution <- myVars$DETECTORS$resolution
  # center <- TRUE
  
  detector.r1 <- raster(extent(data), resolution = resolution, crs = proj4string(data))
  maindetector.r <- rasterize(data, detector.r1)
  
  ### ==== CREATE POLYGONS FROM RASTER  ====
  ## Main detectors 
  temp.r <- raster(maindetector.r)
  temp.r[] <- 1:length(temp.r)
  maindetector.poly <- rasterToPolygons(temp.r, dissolve = TRUE)
  
  ### ==== OBTAIN SPATIALPOINTS FROM DETECTORS ====
  ## Main detectors 
  main.detector.xy <- xyFromCell(maindetector.r, 1:ncell(maindetector.r))
  main.detector.sp <- SpatialPointsDataFrame( data.frame(main.detector.xy[,c("x","y")]),
                                              data=data.frame(main.detector.xy),
                                              proj4string=CRS(projection(data)))
  names(main.detector.sp@data) <- c("main.cell.x","main.cell.y")
  main.detector.sp@data$main.cell.id <- 1:length(main.detector.sp)
  
  myDetectors <- list( detector.sp = main.detector.sp, 
                       grid.poly = maindetector.poly)
  return(myDetectors)
}

###############################################################
##################### UTMToGrid ###############################
###############################################################
#' @title Function to convert UTM coordinates to GRID coordinates for both habitat and detectors
#'
#' @description
#' \code{UTMToGrid} returns a list object with \code{} 
#' 
#' @param data.sp  A \code{SpatialPointsDataFrame} object with detectors' locations.
#' @param grid.sp A \code{SpatialPointsDataFrame} object with habitat grid cell locations.

#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' scaled<-UTMToGrid(grid.sp=habitat.sp,detector.sp=detector.sp)

UTMToGrid <- function(data.sp = NULL,
                      grid.sp = NULL){
  # PREPARE THE DATA
  grid.xy <- as.array(coordinates(grid.sp))
  dimnames(grid.xy) <- list(1:length(grid.sp), c("x","y"))
  data.xy <- as.array(coordinates(data.sp))
  dimnames(data.xy) <- list(1:length(data.sp), c("x","y"))
  # CALCULATE THE RESOLUTION
  resolution <- min(diff(unique(sort(grid.xy[ ,"x"]))))#---assumes square grid cells and utm projection (units in meters or km; not latlong!)
  ## obtain x and y min
  start0.y <- max(grid.xy[ ,"y"]) + resolution/2 #---because we are moving from top to bottom
  start0.x <- min(grid.xy[ ,"x"]) - resolution/2 #---because we are moving from left to right
  ##---- Re-projecting the grid cell centers
  grid.scaled.xy <- grid.xy
  grid.scaled.xy[ ,"y"] <- (start0.y - grid.xy[ ,"y"])/resolution
  grid.scaled.xy[ ,"x"] <- (grid.xy[ ,"x"] - start0.x)/resolution
  ##---- REPROJECTING THE DATA
  data.scaled.xy <- data.xy
  data.scaled.xy[ ,"y"] <- (start0.y - data.xy[ ,"y"])/resolution
  data.scaled.xy[ ,"x"] <- (data.xy[ ,"x"] - start0.x)/resolution 
  out <- list(grid.scaled.xy = grid.scaled.xy,
              grid.xy = grid.xy,
              data.scaled.xy = data.scaled.xy,
              data.xy = data.xy)
  return(out)
}

###############################################################
##################### unSparse_yCombined ###############################
###############################################################

#' Unsparse a Sparse Matrix Preparation
#'
#' R utility function to turn a sparse matrix representation into a two -dimensional detection array. 
#' 
#' Simulation from a nimble model where the \code{\link{dbinomLocal_normal}} and \code{\link{dpois_sparseLocalSCR}}
#' functions are used, returns detection array in sparse matrix format 'yCombined'. We use 
#' \code{unSparse_yCombined} function to  get back the original 'y' matrix that was later 
#' transformed to sparse format to build the nimble model (using \code{getSparseY} function),  
#' advance of model building to create a sparse matrix representation of the observation data. 
#' It creates and returns a matrix:
#' 
#'
#' @param yCombined A two dimensional observation data array with dimension n.individuals x lengthYCombined in sparse format obtained from \code{getSparseY} function.
#' @param ntraps The number of traps in the detector array. 
#'
#' 
#' @return a matrix object represnting the observation data:
#' 
#' \itemize{
#' \item \emph{y} an array of dimensions n.individuals x ntraps  which contains the SCR observations of each individual at each detector.
#' at the traps it was detected at.
#' }
#'
#' @author Soumen Dey
#'
#' @examples
#' y.full <- matrix(rbinom(5000, 5, 0.02), ncol = 100)
#' yCombined <- nimbleSCR::getSparseY(y.full)$yCombined[, , 1]
#' y <- unSparse_yCombined(yCombined, ntraps = 100)
#' all.equal(y, y.full) # TRUE
#' 
#' @export
#' 
#' 
unSparse_yCombined <- function(yCombined, ntraps){
  
  M <- dim(yCombined)[1]
  MAX <- (dim(yCombined)[2] - 1)/2
  ysparse <- yCombined[,2:(MAX+1)]
  DetIndices <- yCombined[,(MAX+2):(2*MAX+1)]
  DetNums <- yCombined[,1]
  y <- do.call(rbind, lapply(1:M, function(ii){
    ## RECREATE THE FULL LENGTH VECTOR OF DETECTIONS
    out <- rep(0, length = ntraps)
    if(DetNums[ii] > 0){
      for(j in 1:DetNums[ii]){
        out[DetIndices[ii,j]] <- ysparse[ii,j]
      }#j
    }#if
    return(out)
  }))
  return(y)
}

# ### NEED THIS DURING DETECTION SIMULATION
# inv.logit <- nimbleFunction(run = function(inValues = double(1)) {
#   a <- 1.0 / (1.0 + exp(-inValues))
#   returnType(double(1))
#   return(a)
# })

###############################################################
##################### e2dist ###############################
###############################################################

e2dist <- function(x,y){ # m x n # 'x' is matrix mx2, 'y' is matrix nx2
  if(!is.matrix(x)){
    d2vec = (x[1] - y[,1])^2 + (x[2] - y[,2])^2 
    return(d2vec)
  }
  if(is.matrix(x)){
    indices <- rep(1:nrow(y),each=nrow(x))
    d2vec = (x[,1] - y[indices,1])^2 + (x[,2] - y[indices,2])^2 
    return(matrix(d2vec, nrow = nrow(x), ncol = nrow(y), byrow = F))
  }
}

###############################################################
##################### fun.trapsInZones ###############################
###############################################################

fun.trapsInZones <- function(n.detectors, ntrapsPerZone){
  # ntrapsPerZone <- n.detectors/nZones
  nZones <- n.detectors/ntrapsPerZone
  #=== Starting point of each zone in one row of zones
  seq1 <-seq(1,sqrt(n.detectors),by = sqrt(ntrapsPerZone))
  #=== Starting point minus 1 for each row of zones
  seq2 <- seq(0, n.detectors, by = sqrt(n.detectors)*sqrt(ntrapsPerZone))
  seq2 <- seq2[-length(seq2)] # Not counting last.detector
  #=== Starting point of each zone in the detector array
  seq3 <- c()
  for(j in 1:length(seq2)){seq3 <- c(seq3, seq2[j]+seq1)}
  #=== Detectors for one zone
  aseq <- c()
  for(j2 in 1:sqrt(ntrapsPerZone)){
    aseq <- c(aseq, (j2-1)*sqrt(n.detectors)+ c(0:(sqrt(ntrapsPerZone)-1)))
  } 
  #=== Detectors in each zone
  trapsInZones <- matrix(NA, ntrapsPerZone, nZones)
  for(j in 1:length(seq3)){ trapsInZones[,j] <- seq3[j] + aseq}
  return(trapsInZones)
}

###############################################################
##################### MakeAugmentation ########################
###############################################################
#' @title MakeAugmentation function
#'
#' @description \code{MakeAugmentation} increases the dimensions of an object along
#'  the individual and/or the year dimension. It returns a \code{Vector} or \code{Matrix} 
#'  object with the expanded object.
#'
#' @param y A \code{Vector} or \code{Matrix} object containing the individual detection histories.
#' @param aug.factor A \code{numeric} object defining the augmentation factor to be used.
#' @param replace.value A \code{numeric} object defining the value to be repeated for augmented individuals.
#' 
#' @return A \code{Vector} or \code{Matrix} object containing the augmented y.

MakeAugmentation <- function( y,
                              aug.factor= NULL,
                              aug.years = NULL,
                              replace.value = NA){
  ## Vector Data augmentation
  if(is.vector(y)){
    if(is.null(names(y))){ 
      names(y) <- 1:length(y) 
    }
    if(!is.null(aug.factor)){
      y.aug <- c(y, rep(replace.value, round(length(y) * aug.factor)))
      names(y.aug) <- c(names(y), rep("Augmented", round(length(y) * aug.factor)))
      y <- y.aug
    }
  }
  
  ## Matrix augmentation
  if(is.matrix(y)){
    if(is.null(dimnames(y))){
      dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2])
    }
    if(!is.null(aug.factor)){
      n.tot <- round(dim(y)[1]*(1 + aug.factor))
      y.aug <- matrix(replace.value, n.tot, dim(y)[2])
      y.aug[1:dim(y)[1], ] <- y
      dimnames(y.aug) <- list(c( dimnames(y)[[1]], rep("Augmented", n.tot - dim(y)[1])),
                              dimnames(y)[[2]])
      y <- y.aug
    }
    
  }  
  
  return (y)
}


###############################################################
##################### MakeInitXY ##############################
###############################################################

#' @title Function to set initial values of XY coordinates of ACS
#'
#' @description
#' \code{InitXY} returns a matrix object with with the x coordinates ([,1]) and y coordinates  ([,2]). it returns the location of a detection for an detected individual and a random location for augmented individuals. 
#' 
#' @param y \code{matrix}  with individual detections. row=Ids, col= detectors. 
#' @param habitat.mx \code{matrix} of corrdinates from the buffered habitat grid
#' @param detector.xy \code{matrix} with coordinates of detectors scaled
#' @param IDCells.mx \code{matrix} the ID of the cells of the habitat
#' @param grid.xy \code{matrix} with the scaled coordinates of the habitat 
#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' MakeInitXY(   y = data$y
#'              , habitat.mx = data$habitat.mx
#'              , detector.xy = data$detector.xy
#'              , grid.xy = grid.xy)

MakeInitXY <- function(    y = y
                           , habitat.mx = habitat.mx
                           , detector.xy = detector.xy
                           , IDCells.mx = IDCells.mx
                           , grid.xy = grid.xy){
  # , xy.bounds = NULL ){
  
  # ---- STEP 1: GET THE OBJECT READY TO STORE THE DATA ----- 
  # Make it general to work with the time dimention 
  n.years <- 1
  # n.years <- ifelse(length(dim(y))>=3, dim(y)[3], 1)
  n.individuals <- dim(y)[1]
  n.detectors <- dim(detector.xy)[1]
  
  # if(length(dim(y)) == 2){y <- array(y, c(n.individuals, n.detectors, n.years))}
  # if(length(dim(detector.xy)) == 2){detector.xy <- array(detector.xy, c(n.detectors, 2, n.years))}
  # n.year <- dim(y)[3]
  
  
  # STORE THE DIMENSION OF THE WINDOW 
  # IF NO FRAGMENTATION IS USED 
  dim.mx <- array(0, c(n.individuals, 2, 2))#, n.year))
  # for(t in 1:n.year){
  dim.mx[,1,1] <- 1 # min x
  dim.mx[,2,1] <- 1 # min y
  dim.mx[,1,2] <- dim(habitat.mx)[2]# max x
  dim.mx[,2,2] <- dim(habitat.mx)[1]# max y 
  
  # if(!is.null(xy.bounds)){
  #     dim.mx <- xy.bounds
  #     dim.mx[,1,1] <-  floor(dim.mx[,1,1]) + 1 # Add one to match the habitat 
  #     dim.mx[,2,1] <-  floor(dim.mx[,2,1]) + 1 # Add one to match the habitat
  # }
  # }
  
  # IF  FRAGMENTATION IS USED 
  
  
  
  # empty <- array(numeric(),c(n.individuals, 2, n.year))
  empty <- array(numeric(),c(n.individuals, 2))
  ids <- lapply(c(1:n.individuals), function(x) x)
  
  # ---- STEP 2: FIND SXY FOR ALL INDIVIDUALS  ----- 
  
  # for(t in 1:n.year){
  listt <- apply(y, 1, function(x) which(x>0, arr.ind = T))# obtain detectors with detections for each ids 
  
  if(length(listt)==0){listt <- list()
  for(i in 1:n.individuals){
    listt[[i]] <- integer()
  }
  }
  
  # empty[,,t] <- do.call(rbind, lapply(ids, function(i) {
  empty <- do.call(rbind, lapply(ids, function(i) {
    
    #print(i)
    x <-  listt[[i]]
    # if only one detections 
    if(length(x)==1){ # If only one detection, use that detection as a starting value
      detector.xy[x,]
    }else{
      
      # if several detections    
      if(length(x)>1){# if more than 1 detection, use average value
        mean.det <- colMeans(detector.xy[x,])
        # if doesnt end up in habitat, use a random coordinate within habitat 
        if(habitat.mx[floor(mean.det)[2]+1,floor(mean.det)[1]+1]==1){
          mean.det 
        }else{
          min.x <- dim.mx[i,1,1]
          max.x <- dim.mx[i,1,2]
          min.y <- dim.mx[i,2,1]
          max.y <- dim.mx[i,2,2]
          
          window.mx.ind <- habitat.mx[c(min.y:max.y), c(min.x:max.x)]
          window.ID.ind <- IDCells.mx[c(min.y:max.y), c(min.x:max.x)]
          
          id.cell <- window.ID.ind[window.mx.ind==1]
          sxy.grid <- matrix(grid.xy[id.cell,], nrow=length(id.cell),ncol=2)
          
          matrix(sxy.grid[floor(runif(1, 1, length(id.cell))),], ncol=2, nrow=1)                     }
      }else{
        
        # if no detections (augmented IDS)
        if(length(x)==0){
          
          min.x <- dim.mx[i,1,1]
          max.x <- dim.mx[i,1,2]
          min.y <- dim.mx[i,2,1]
          max.y <- dim.mx[i,2,2]
          
          window.mx.ind <- habitat.mx[c(min.y:max.y), c(min.x:max.x)]
          window.ID.ind <- IDCells.mx[c(min.y:max.y), c(min.x:max.x)]
          
          id.cell <- window.ID.ind[window.mx.ind==1]
          sxy.grid <- matrix(grid.xy[id.cell,], nrow=length(id.cell),ncol=2)
          
          matrix(sxy.grid[floor(runif(1, 1, length(id.cell))),], ncol=2, nrow=1)
        }
      }}      
  })) 
  
  # }
  return(empty)
  
}
###############################################################
##################### dbern_vector_noIndicator ##############################
###############################################################

#' @title Function to create a NIMBLE custom distribution with internalized detection
#'  probabilities calculation for faster SCR model runs.
#'
#' @description
#' \code{dbern_vector_noIndicator} returns the likelihood of a given binary vector x[i,1:length.x] 
#' 
#' @param x \code{Vector} containing observation/non-observations 
#' @param prob \code{Numeric} variable denoting the success probability (i.e., 1).
#' @param length.x A \code{Numeric} denoting the length of x. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)
#'
#' @examples
#' x[i,1:n.detectors] ~ dbern_vector_noIndicator(p0 , prob, length.x)

#### 1.Density function ####
dbern_vector_noIndicator <- nimbleFunction(run = function( x = double(1)
                                                           , prob = double(0)
                                                           , length.x = double(0) #length of x
                                                           , log = integer(0, default = 0)){
  # Return type declaration
  returnType(double(0))
  
  
  ## Calculate the likelihood of the detection observations
  logProb <- 0
  for(j in 1:length.x){
    # Calculate the log-likelihood of each observation
    # logProb <- logProb + dbern(x[j], prob = prob, log = TRUE)
    if(x[j] == 0){
      logProb <- logProb + log(1.0 - prob)
    }else{ logProb <- logProb + log(prob)}
  }#j
  
  # Output
  if(log)return(logProb)
  return(exp(logProb))
  
})

#### 2.Sampling function ####
rbern_vector_noIndicator <- nimbleFunction(run = function( n = integer(0)
                                                           , prob = double(0)
                                                           , length.x = double(0)
                                                           # , indicator = double(0, default = 1.0)
){
  # Return type declaration
  returnType(double(1))
  
  # Check input dimensions
  if(n!=1){print("rbern_vector only allows n = 1; using n = 1")}
  
  # Draw from a Bernoulli distribution with the calculated probability
  detectOut <- rbinom(length.x, 1, prob)
  
  # Output
  return(detectOut)
})

#### 3.Registration ####
registerDistributions(list(
  dbern_vector_noIndicator = list(
    BUGSdist = "dbern_vector_noIndicator(prob, length.x)",
    types = c( "value = double(1)", "prob = double(0)", "length.x = double(0)"),
    pqAvail = FALSE)))

####################################################################
### binary Gibbs sampler for vector parameter ###########################################
####################################################################

#' @rdname samplers
#' @export

sampler_binaryFEWatAtime <- nimbleFunction(
  name = 'sampler_binaryFEWatAtime',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    
    numIndicesToChange <- extractControlElement(control, "numIndicesToChange", 2)   ## default value is 2
    
    length.target <- length(targetAsScalar)
    prop.value <- numeric(length.target)
    probVec <- rep(1/length.target, length.target)
    thisTwoIndices <- numeric(2)
    ## checks
    # if(length(targetAsScalar) > 1)  stop('cannot use binary sampler on more than one target node')
    #####if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
  },
  run = function() {
    currentLogProb <- model$getLogProb(calcNodes)
    prop.value <<- model[[target]]
    for(i in 1:numIndicesToChange) {
      changeIndex <<- rcat(n = 1, prob = probVec)
      prop.value[changeIndex] <<- 1 - prop.value[changeIndex]
    }
    model[[target]] <<- prop.value
    # model[[target]] <<- 1 - model[[target]]
    otherLogProbPrior <- model$calculate(target)
    if(otherLogProbPrior == -Inf) {
      otherLogProb <- otherLogProbPrior
    } else {
      otherLogProb <- otherLogProbPrior + model$calculate(calcNodesNoSelf)
    }
    acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1) # = otherProb / (currentProb+otherProb)
    ##  The runif(1,0,1) < acceptanceProb in below is not a MH step, 
    ##  this is just alternative computation, Bernoulli full conditional
    jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb) 
    if(jump) {
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    }
  },
  methods = list(
    reset = function() { }#,
    # indicesToggle = function() {
    #   this.indices <- sample(1:length(target), 2, replace=FALSE)
    #   returnType(double(1))
    #   return(this.indices)
    # }
  )
)


########################################################################
##-- BASIC SINGLE-SEASON SCR
########################################################################
modelCode_SCR <- nimbleCode({
  
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  ##---- DETECTION PROCESS
  sigma ~ dunif(0, 50)
  logit.p0 ~ dnorm(mean = 0, sd = 2)
  p0 <- ilogit(logit.p0)
  p0Vec[1:n.detectors] <- p0
  for (i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
    z[i] ~ dbern(psi)
    
    yBinded[i, 1:dimYbinded] ~ dbinomLocal_normal(
      size = trials[1:n.detectors],
      p0 = p0,
      s = sxy[i,1:2],
      sigma = sigma,
      trapCoords =  detector.xy[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.cells,1:maxNBDets],
      localTrapsNum = numLocalTraps[1:n.cells],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatGridDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      lengthYCombined = dimYbinded)
    
  }#i
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})


########################################################################
##-- SCR-RE - RANDOM NORMAL EFFECT ON DETECTORS
########################################################################

modelCode_RE <- nimbleCode({
  
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0, 50)
  sigma.w ~ dgamma(1, 1)
  
  logit.eta ~ dnorm(mean = 0, sd = 2)
  eta <- ilogit(logit.eta)
  mu[1:n.detectors] <- rep(0, n.detectors)
  chol.v.mat[1:n.detectors, 1:n.detectors] <- sigma.w*diag(n.detectors) 
  W[1:n.detectors] ~ dmnorm(mean = mu[1:n.detectors],
                                   cholesky = chol.v.mat[1:n.detectors, 1:n.detectors],
                                   prec_param = 0)
  logit(p0Vec[1:n.detectors]) <- logit.eta + W[1:n.detectors]
  
  for (i in 1:n.individuals){
    
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
    z[i] ~ dbern(psi)
    
    yBinded[i, 1:dimYbinded] ~ dbinomLocal_normal(
      size = trials[1:n.detectors],
      p0Traps = p0Vec[1:n.detectors],
      s = sxy[i,1:2],
      sigma = sigma,
      trapCoords =  detector.xy[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.cells,1:maxNBDets],
      localTrapsNum = numLocalTraps[1:n.cells],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatGridDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      lengthYCombined = dimYbinded)
    
  }#i
  
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})
########################################################################
##-- SCR-RE.agg - RANDOM NORMAL EFFECT ON DETECTORS - GROUPED INTO CLUSTERS
########################################################################
modelCode_RE.agg <- nimbleCode({
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0, 50)
  sigma.w ~ dgamma(1, 1)
  logit.eta ~ dnorm(mean = 0, sd = 2)
  eta <- ilogit(logit.eta)
  mu[1:nZones] <- rep(0, nZones)
  chol.v.mat[1:nZones, 1:nZones] <- sigma.w*diag(nZones) 
  Wzone[1:nZones] ~ dmnorm(mean = mu[1:nZones],
                                cholesky = chol.v.mat[1:nZones, 1:nZones],
                                prec_param = 0)
  W[1:n.detectors] <- B[1:n.detectors,1:nZones] %*% Wzone[1:nZones]
  logit(p0Vec[1:n.detectors]) <- logit.eta + W[1:n.detectors]
  
  
  for (i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
    z[i] ~ dbern(psi)
    
    yBinded[i, 1:dimYbinded] ~ dbinomLocal_normal(
      size = trials[1:n.detectors],
      p0Traps = p0Vec[1:n.detectors],
      s = sxy[i,1:2],
      sigma = sigma,
      trapCoords =  detector.xy[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.cells,1:maxNBDets],
      localTrapsNum = numLocalTraps[1:n.cells],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatGridDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      lengthYCombined = dimYbinded)
  }#i
  
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})
########################################################################
##-- SCR-SARE - Spatially autocorrelated random effect 
########################################################################
modelCode_SARE <- nimbleCode({
  logit.eta ~ dnorm(mean = 0, sd = 2)
  eta <- ilogit(logit.eta)
  sigma ~ dunif(0,50)
  log.phi ~ dnorm(mean = 0, sd = 5)
  phi <- exp(log.phi)
  psi ~ dunif(0,1)
  v.mat[1:n.detectors, 1:n.detectors] <- exp(- phi * distance[1:n.detectors,1:n.detectors])
  chol.v.mat[1:n.detectors, 1:n.detectors] <- chol(v.mat[1:n.detectors,1:n.detectors])
  
  mu[1:n.detectors] <- logit.eta
  W[1:n.detectors] ~ dmnorm(
    mean = mu[1:n.detectors],
    cholesky = chol.v.mat[1:n.detectors,1:n.detectors],
    prec_param = 0)
  logit(p0Vec[1:n.detectors]) <- W[1:n.detectors]
  
  
  for (i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
    z[i] ~ dbern(psi)
    yBinded[i, 1:dimYbinded] ~ dbinomLocal_normal(
      size = trials[1:n.detectors],
      p0Traps = p0Vec[1:n.detectors],
      s = sxy[i,1:2],
      sigma = sigma,
      trapCoords =  detector.xy[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.cells,1:maxNBDets],
      localTrapsNum = numLocalTraps[1:n.cells],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatGridDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      lengthYCombined = dimYbinded)
  }#i
  N <- sum(z[1:n.individuals])
})


########################################################################
##-- SCR-SARE.agg - Spatially autocorrelated random effect - GROUPED INTO CLUSTERS
########################################################################
modelCode_SARE.agg <- nimbleCode({
  
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0, 50)
  logit.eta ~ dnorm(mean = 0, sd = 2)
  eta <- ilogit(logit.eta)
  
  log.phi ~ dnorm(mean = 0, sd = 5)
  phi <- exp(log.phi)
  
  v.mat[1:nZones, 1:nZones] <- exp(- phi * distance.zone[1:nZones,1:nZones])
  chol.v.mat[1:nZones, 1:nZones] <- chol(v.mat[1:nZones,1:nZones])
  
  mu[1:nZones] <- logit.eta
  Wzone[1:nZones] ~ dmnorm(
    mean = mu[1:nZones],
    cholesky = chol.v.mat[1:nZones,1:nZones],
    prec_param = 0)
  
  
  W[1:n.detectors] <- B[1:n.detectors,1:nZones] %*% Wzone[1:nZones]
  logit(p0Vec[1:n.detectors]) <- W[1:n.detectors]
  
  for (i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
    z[i] ~ dbern(psi)
    
    yBinded[i, 1:dimYbinded] ~ dbinomLocal_normal(
      size = trials[1:n.detectors],
      p0Traps = p0Vec[1:n.detectors],
      s = sxy[i,1:2],
      sigma = sigma,
      trapCoords =  detector.xy[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.cells,1:maxNBDets],
      localTrapsNum = numLocalTraps[1:n.cells],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatGridDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      lengthYCombined = dimYbinded)
  }#i
  
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})

########################################################################
##-- SCR_FM - 2-Group Finite Mixture Capture-Recapture Model
########################################################################

modelCode_FM <- nimbleCode({
  
  eta[1] ~ dunif(0, 1)
  eta[2] ~ dunif(0, 1)
  one ~ dconstraint(eta[1] <= eta[2])
  pi ~ dbeta(1,1) 
  u[1:n.detectors] ~ dbern_vector_noIndicator(prob = pi, length.x = n.detectors) 
  p0Vec[1:n.detectors] <- (1-u[1:n.detectors])*eta[1] +  u[1:n.detectors]*eta[2]
  sigma ~ dunif(0,50)
  psi ~ dunif(0,1)
  
  for (i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
    z[i] ~ dbern(psi)
    yBinded[i, 1:dimYbinded] ~ dbinomLocal_normal(
      size = trials[1:n.detectors],
      p0Traps = p0Vec[1:n.detectors],
      s = sxy[i,1:2],
      sigma = sigma,
      trapCoords =  detector.xy[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.cells,1:maxNBDets],
      localTrapsNum = numLocalTraps[1:n.cells],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatGridDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      lengthYCombined = dimYbinded)
  }#i
  N <- sum(z[1:n.individuals])
  
})


########################################################################
##-- SCR_FM.agg - 2-Group Finite Mixture Capture-Recapture Model - GROUPED INTO CLUSTERS
########################################################################
modelCode_FM.agg <- nimbleCode({
  
  eta[1] ~ dunif(0, 1)
  eta[2] ~ dunif(0, 1)
  one ~ dconstraint(eta[1] <= eta[2])
  pi ~ dunif(0, 1) 
  for (this.zone in 1:nZones){
    u[this.zone] ~ dbern(pi) 
  }
  p0Zone[1:nZones] <- (1-u[1:nZones])*eta[1] +  u[1:nZones]*eta[2]
  p0Vec[1:n.detectors] <- B[1:n.detectors,1:nZones] %*% p0Zone[1:nZones]
  sigma ~ dunif(0,50)
  psi ~ dunif(0,1)
  
  for (i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
    z[i] ~ dbern(psi)
    yBinded[i, 1:dimYbinded] ~ dbinomLocal_normal(
      size = trials[1:n.detectors],
      p0Traps = p0Vec[1:n.detectors],
      s = sxy[i,1:2],
      sigma = sigma,
      trapCoords =  detector.xy[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.cells,1:maxNBDets],
      localTrapsNum = numLocalTraps[1:n.cells],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatGridDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      lengthYCombined = dimYbinded)
  }#i
  N <- sum(z[1:n.individuals])
  
})


########################################################################
##-- FE - KNOWN TRUE COVARIATE ON DETECTORS
########################################################################
modelCode_FE <- nimbleCode({
  logit.eta ~ dnorm(mean = 0, sd = 2)
  eta <- ilogit(logit.eta)
  sigma ~ dunif(0,50)
  log.phi ~ dnorm(mean = 0, sd = 5)
  phi <- exp(log.phi)
  psi ~ dunif(0,1)
  v.mat[1:n.detectors, 1:n.detectors] <- exp(- phi * distance[1:n.detectors,1:n.detectors])
  chol.v.mat[1:n.detectors, 1:n.detectors] <- chol(v.mat[1:n.detectors,1:n.detectors])
  
  mu[1:n.detectors] <- logit.eta
  W[1:n.detectors] ~ dmnorm(
    mean = mu[1:n.detectors],
    cholesky = chol.v.mat[1:n.detectors,1:n.detectors],
    prec_param = 0)
  logit(p0Vec[1:n.detectors]) <- W[1:n.detectors]
  
  for (i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
    z[i] ~ dbern(psi)
    yBinded[i, 1:dimYbinded] ~ dbinomLocal_normal(
      size = trials[1:n.detectors],
      p0Traps = p0Vec[1:n.detectors],
      s = sxy[i,1:2],
      sigma = sigma,
      trapCoords =  detector.xy[1:n.detectors,1:2],
      localTrapsIndices = localTrapsIndices[1:n.cells,1:maxNBDets],
      localTrapsNum = numLocalTraps[1:n.cells],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatGridDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],# 1,
      lengthYCombined = dimYbinded)
  }#i
  N <- sum(z[1:n.individuals])
  
})

