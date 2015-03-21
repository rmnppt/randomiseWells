# Randomise treatments on a 96 well plate

# The first function, randomiseWells() takes some information about the plate and samples 
# and then does the randomisation. 
randomiseWells <- function(nRows, nCols, treatments, sampleSizes = NULL, emptyOuter = 0){
  
  # outer wells empty?
  if(emptyOuter > 0){nRows <- nRows-2*emptyOuter; nCols <- nCols-2*emptyOuter}
  nWells <- nRows * nCols
  
  # no sample sizes given?
  if(is.null(sampleSizes)){sampleSizes <- rep(floor(nWells / length(treatments)), length(treatments))}
  
  # account for empty space
  emptySpace <- nWells - sum(rep(1, length(treatments)) * sampleSizes)
  if(emptySpace > 0){
    treatments <- c(treatments, NA)
    sampleSizes <- c(sampleSizes, emptySpace)
  }
  
  # calculate the pool
  pool <- numeric(0)
  for(i in 1:length(treatments)){
    pool <- c(pool, rep(treatments[i], sampleSizes[i]))
  }
  
  # pick random treatments for wells
  randomPool <- sample(pool, replace = F)
  plate <- matrix(randomPool, nRows, nCols)
  if(emptyOuter > 0){for(i in 1:emptyOuter){
    plate <- rbind(NA, cbind(NA, plate, NA), NA)
  }}
  
  return(plate)
}


# The second function randomNeighbours() takes a matrix object (typically returned by randomiseWells())
# and records all neighbouring elements for each element of the matrix. Then it peforms a CHI squared
# test to see if neighbours are distributed equally. This will not work when a sufficiently large number
# of treatments is spread over a sufficiently small plate.
randomNeighbours <- function(plate){
  
  nRows <- dim(plate)[1]
  nCols <- dim(plate)[2]
  
  # first define the identities and neighbours 
  identities <- c(plate)
  paddedPlate <- rbind(NA, cbind(NA, plate, NA), NA)
  neigh <- rbind(N  = as.vector(paddedPlate[1:nRows, 2:(nCols + 1)]),
                 NE = as.vector(paddedPlate[1:nRows, 3:(nCols + 2)]),
                 E  = as.vector(paddedPlate[2:(nRows + 1), 3:(nCols + 2)]),
                 SE = as.vector(paddedPlate[3:(nRows + 2), 3:(nCols + 2)]),
                 S  = as.vector(paddedPlate[3:(nRows + 2), 2:(nCols + 1)]),
                 SW = as.vector(paddedPlate[3:(nRows + 2), 1:nCols]),
                 W  = as.vector(paddedPlate[2:(nRows + 1), 1:nCols]),
                 NW = as.vector(paddedPlate[1:nRows, 1:nCols]))
  
  # now count up the number of each treatment surrounding each treatment
  treatments <- unique(identities)
  treatments <- treatments[which(!is.na(treatments))]

  sampleSizes <- numeric(0)
  for(i in 1:length(treatments)){
    sampleSizes[i] <- length(which(identities == treatments[i]))
  }
    
  neighData <- matrix(NA, length(treatments), length(treatments), 
                      dimnames = list(treatments,treatments))
  expected <- neighData
  totalNeighbours <- length(which(!is.na(neigh)))
  for(i in 1:length(treatments)){
    these <- c(neigh[,which(identities == treatments[i])])
    for(j in 1:length(treatments)){
      neighData[j,i] <- length(which(these == treatments[j]))
    }
    expected[,i] <- totalNeighbours*(sampleSizes/nWells)*(sampleSizes[i]/nWells)
  }
  
  # perfor a chisq.test to test whether sufficiently 'random'
  these <- which(upper.tri(expected, T))
  return(chisq.test(neighData[these], p = expected[these], rescale.p = T))
  
}


# The third function edgeTest() takes a matrix object (typically produced with randomiseWells()),
# and measures the distribution of treatments at the edge of the plate.
edgeTest <- function(plate){
  
  nRows <- dim(plate)[1]
  nCols <- dim(plate)[2]
  # identify treatments
  treatments <- unique(c(plate))
  treatments <- treatments[which(!is.na(treatments))]
  
  # count sample sizes
  sampleSizes <- numeric(0)
  for(i in 1:length(treatments)){
    sampleSizes[i] <- length(which(identities == treatments[i]))
  }
  
  # if there are empty outer wells
  if(all(is.na(plate[1,]))){ 
    # remove them
    plate <- plate[-c(1, nRows), -c(1, nCols)]
    # and update size
    nRows <- dim(plate)[1]
    nCols <- dim(plate)[2]
  }
  
  # identify outer wells
  outerWells <- c(
    plate[1,], plate[nRows,], plate[-c(1, nRows), 1], plate[-c(1, nRows), nCols]
  )
  
  # count em up
  outerCount <- numeric(0)
  for(i in 1:length(treatments)){
    outerCount[i] <- length(which(outerWells == treatments[i]))
  }

  return(chisq.test(outerCount, p = sampleSizes, rescale.p = T))
}




