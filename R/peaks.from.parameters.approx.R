
peaksFromParametersApprox <- function(aCVec, a, b, c, d, z, phi, q0, nrRRoots, nrCRoots, stopOption, nrPeaks, coverage, abundantEstim, approx, approxStart, approxParam){
##used in calculatePeaks (defined below) and calculateDifferential (calculate.differential.R)
  
 
  if (is.null(nrPeaks)){      
      aC <- getAC(aCVec)      
      nrPeaks <- calculateNrPeaks(aC)      
  }
  
  if (!is.na(pmatch(stopOption, "nrPeaks"))) 
        stopOption <- "nrPeaks"
#   print(stopOption)
  STOPS <- c("nrPeaks", "coverage", "abundantEstim")
  stopIdx <- pmatch(stopOption, STOPS)
#   print(stopOption)
  
  if (is.na(stopIdx)) 
      stop("invalid stop option")
  if (stopIdx == -1) 
      stop("ambiguous stop option")

  if (is.na(nrPeaks))
      stop("specify maximal number of calculated peaks")


  if ((stopOption == "coverage"))
      if (is.null(coverage))
	stop("this stop option is not allowed for approx=TRUE")

  if (stopOption == "abundantEstim")
      if (is.null(abundantEstim))
	stop("this stop option is not allowed for approx=TRUE")

  if (is.null(approxParam))
    approxParam = nrPeaks		
	#stop("approxParam can not be NULL for approx=TRUE")

  v.r <- rep(aCVec, nrRRoots)
  v.c <- rep(aCVec, nrCRoots)
  
  A.r.tmp <- -v.r 
  A.c.tmp <- - 2 * v.c


  calculate.A.r <- function(k){
  l <- length(v.r)
  if (l > 0){
    sum(A.r.tmp)
    }else{
      0
    }  
  }

  calculate.A.c <- function(k){ 
  l <- length(v.c)
  if (l > 0){
    sum(A.c.tmp * cos(k * phi)) 
  }else{
    0
  }  
  }

  calculate.A <- function(k){
    calculate.A.r(k) + calculate.A.c(k)
  }

  ### start of the body of peaks.from.parameters
  ### declarations


 start <- ifelse(approxStart <= nrPeaks, approxStart, 1) 

#  start <- max(approxStart, nrPeaks-approxStart+1) #nrPeaks - approxParam 
#print("start:")
#print(start)
#print(nrPeaks)
#print(nrPeaks-start+1)
  i <- 1


  peaks <- numeric(nrPeaks-start+1) #TODO: #cyclic buffer for peaks table
  A <- numeric(approxParam)

  if (start == 1){	
   peaks[i] <- exp(sum(aCVec*log(q0)))
  }
  else{
   peaks[i] <- 1	
  }	 	
	
	
 for (i in 2:(approxParam+1)){
    A.r.tmp <- A.r.tmp / b
    A.c.tmp <- A.c.tmp / sqrt(c^2 + d^2)
    A[i-1] <- calculate.A(i-1)    ##przerzucic liczenie tej tablicy do osobnej - stalej - petli
}

  ### parameters initialization

	
  
  ### main loop for calculating peaks
  i <- 2
  length <- 1

  while (i <= nrPeaks-start+1){     
#   for (i in 2:nr.peaks){

    low.idx <- max(1, i-approxParam)
    hi.idx <- min(i-1, approxParam)
    peaks[i] <- sum(peaks[low.idx:(i-1)] * A[hi.idx:1])/(i-1+start-1)#(i-1+start-1)##(i-1) 
    i <- checkOption(peaks, i, stopOption, nrPeaks, coverage, abundantEstim)    
    i <- i + 1
    length <- length + 1   
  } 
  peaks[1:length]
  ### result stored in 'peaks'
}


