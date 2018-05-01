#' @useDynLib quasar
#' @importFrom Rcpp evalCpp


#' @export
rConvolve <- function(quasarcl, inputMatrix, sizesVector, filterVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(nrow(inputMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	if(!identical(nrow(inputMatrix), length(filterVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppConvolve(quasarcl, inputMatrix, sizesVector, filterVector))
}



#' @export
rCopyIfNotInf <- function(quasarcl, inputMatrix, filteredSize) 
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppCopyIfNotInf(quasarcl, inputMatrix, adjustSize(filteredSize)))
}



#' @export
rCountIfNotInf <- function(quasarcl, inputMatrix) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppCountIfNotInf(quasarcl, inputMatrix))
}



#' @export
rReglin <- function(quasarcl, xMatrix, yMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(xMatrix), dim(yMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(xMatrix), length(sizesVector)))
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppReglin(quasarcl, xMatrix, yMatrix, sizesVector))
}



#' @export
rReglinYax <- function(quasarcl, xMatrix, yMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(xMatrix), dim(yMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(xMatrix), length(sizesVector)))
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppReglinYax(quasarcl, xMatrix, yMatrix, sizesVector))
}



#' @export
rChisq <- function(quasarcl, fMatrix, yMatrix, errorsMatrix, sizesVector)
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(fMatrix), dim(yMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(fMatrix), dim(errorsMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(fMatrix), length(sizesVector)))
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppChisq(quasarcl, fMatrix, yMatrix, errorsMatrix, sizesVector))
}



#' @export
rTrapz <- function(quasarcl, yMatrix, xMatrix, sizesVector)
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(yMatrix), dim(xMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(yMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppTrapz(quasarcl, yMatrix, xMatrix, sizesVector))
}


#wszystkie ważne piki znaleźć
#' @export
extremefits <- function(quasars, wavelengthsMatrix, continuumsMatrix, sizesVector, outputfit)
{
  vi=c()
  va=c()
  for(q in 1:length(quasars)) {
    if (length(outputfit[[q]]$elementsFits) > 0)
      for (i in 1:length(outputfit[[q]]$elementsFits)) {
        lmin <- wavelengthsMatrix[q,][1]
        lmax <- wavelengthsMatrix[q,][sizesVector[q]]
        lx   <- outputfit[[q]]$elementsFits[[i]]$gaussianParams[2]
        fwhm <- outputfit[[q]]$elementsFits[[i]]$gaussianFWHM
        ew   <- outputfit[[q]]$elementsFits[[i]]$ew
        contin <- continuumsMatrix[q,][1:sizesVector[q]]
        cmin <- contin[1]
        cmax <- contin[sizesVector[q]]
        if (cmin > cmax && lx < lmax && lx > lmin)
          if (fwhm > 0 && fwhm < 10000 && ew < 150) {
            if (length(vi[vi == q]) < 1) {
              vi = c(vi, q)
            }
          } else {
            if (length(va[va == q]) < 1 && length(vi[vi == q]) < 1) {
              va = c(va, q)
            }
          }
      }
  }
  va=c()
  for(q in 1:length(quasars)) {
    if (length(outputfit[[q]]$elementsFits) > 0)
      for (i in 1:length(outputfit[[q]]$elementsFits)) {
        lmin <- wavelengthsMatrix[q,][1]
        lmax <- wavelengthsMatrix[q,][sizesVector[q]]
        lx   <- outputfit[[q]]$elementsFits[[i]]$gaussianParams[2]
        fwhm <- outputfit[[q]]$elementsFits[[i]]$gaussianFWHM
        ew   <- outputfit[[q]]$elementsFits[[i]]$ew
        contin <- continuumsMatrix[q,][1:sizesVector[q]]
        cmin <- contin[1]
        cmax <- contin[sizesVector[q]]
        if (cmin > cmax && lx < lmax && lx > lmin)
          if (fwhm > 0 && fwhm < 10000 && ew < 150) {
            if (length(vi[vi == q]) < 1) {
              vi = c(vi, q)
            }
          } else {
            if (length(va[va == q]) < 1 && length(vi[vi == q]) < 1) {
              va = c(va, q)
            }
          }
      }
  }
  va
}

#' @export
memorysort<-function(){
  library(gdata)
  #ll() # return a dataframe that consists of a variable name as rownames, and class and size (in KB) as columns
  #subset(ll(), KB > 1000) # list of object that have over 1000 KB
  ll()[order(ll()$KB),]
}
