ASTRO_OBJ_SPEC_SIZE <- 4096


#' @export
drawGaussianRawData <- function(element, i, lambda)
{
  if (element$fitParams[[i]][4] > 0) {
    lines(lambda, gaussian(lambda, element$fitParams[[i]][1] , element$fitParams[[i]][2], element$fitParams[[i]][3]), ylim=c(0,element$fitParams[[i]][1]),lwd=2)
  }
}

#' @export
plotGaussian <- function(lambda, a, b, c) 
{
		plot(lambda, gaussian(lambda, a, b, c), ylim=c(0,a))
}

#' @export
lineGaussian <- function(lambda, a, b, c)
{
        line(lambda, gaussian(lambda, a, b, c), ylim=c(0,a))
}

#' @export
drawGaussian <- function(element, lambda)
{
    lines(lambda, gaussian(lambda, element$gaussianParams[1] , element$gaussianParams[2], element$gaussianParams[3]), ylim=c(0,element$gaussianParams[1]),lwd=2)
}

#' @export
drawSpectrum <- function(lambda, spectrum, name)
{
  op <- par(mar=c(5,5,2,1))
  plot(lambda,spectrum[1:length(lambda)],ylim=c(0, max(spectrum)*1.1), type='l', lwd=.8, lty=5, col="gray45", main=name,xlab=expression(paste("wavelength ", "[",ring(A),"]")),ylab=TeX('flux $\\lbrack erg$ $s^{-1}$ $cm^{-2}$ $A^{-1} \\rbrack$'));
}

#' @export
drawSpectrumWithPeaksRawData <- function(picturePath, spectrumsMatrix, wavelengthsMatrix, outputfit, qParams, sizes) {
  for (q in seq(1:length(sizes))) {
    jpeg(file=paste(picturePath, "spectrum_", formatC(q, width=6, flag="0"),".jpg",sep=""),width = 500, height = 500, quality = 55, bg = "white")
    drawSpectrum(spectrumsMatrix[q,][1:sizes[q]], wavelengthsMatrix[q,][1:sizes[q]], paste0("SDSS J",substr(qParams[[q]]$name,5,22)))
    lapply(outputfit[[q]]$elementsFits, drawGaussian, wavelengthsMatrix[q,][1:sizes[q]])
    dev.off();
  }
}

#' @export
drawChosenSpectrumWithPeaksRawData <- function(q, spectrumsMatrix, wavelengthsMatrix, outputfit, qParams, sizes) {
  drawSpectrum(wavelengthsMatrix[q,][1:sizes[q]], spectrumsMatrix[q,][1:sizes[q]], paste0("SDSS J",substr(qParams[[q]]$name,5,22)))
  invisible(lapply(outputfit[[q]]$elementsFits, drawGaussian, wavelengthsMatrix[q,][1:sizes[q]]))
}

#' @export
drawSpectrumWithContinuum <- function(chosen_q, quasars, wavelengthsMatrix, spectrumsMatrix, continuumsMatrix, sizesVector) {
  spectrumsMatrixORIG <- getSpectrumsMatrix(quasars)
  qParams<-getParams(quasars)
  op <- par(mar=c(5,5,2,1))
  plot(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],spectrumsMatrixORIG[chosen_q,1:sizesVector[chosen_q]],type = "l",xlab=expression(paste("wavelength ", "[",ring(A),"]")),ylab=TeX('flux $\\lbrack erg$ $s^{-1}$ $cm^{-2}$ $A^{-1} \\rbrack$'),main=paste0("SDSS J",substr(qParams[[chosen_q]]$name,5,22)),ylim=c(0, max(spectrumsMatrix[chosen_q,1:sizesVector[chosen_q]])*1.1),lwd=.1,col="gray60")
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],spectrumsMatrix[chosen_q,1:sizesVector[chosen_q]])
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],continuumsMatrix[chosen_q,1:sizesVector[chosen_q]],col="gray45")
}

#' @export
drawSpectrumWithoutContinuumIron <- function(chosen_q, quasars, wavelengthsMatrix, spectrumsMatrix, continuumsMatrix, feTemplatesMatrix, sizesVector) {
  spectrumsMatrixORIG <- getSpectrumsMatrix(quasars)
  spectrumsMatrixNoIron <- rMinusMatrix(ptr, spectrumsMatrix, feTemplatesMatrix)
  spectrumsMatrixEmissionLines <- rMinusMatrix(ptr, spectrumsMatrixNoIron, continuumsMatrix)
  qParams<-getParams(quasars)
  op <- par(mar=c(5,5,2,1))
  plot(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],spectrumsMatrix[chosen_q,1:sizesVector[chosen_q]],type = "l",xlab=expression(paste("wavelength ", "[",ring(A),"]")),ylab=TeX('flux $\\lbrack erg$ $s^{-1}$ $cm^{-2}$ $A^{-1} \\rbrack$'),main=paste0("SDSS J",substr(qParams[[chosen_q]]$name,5,22)),ylim=c(0, max(spectrumsMatrix[chosen_q,1:sizesVector[chosen_q]])*1.1),lwd=0.5)
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],continuumsMatrix[chosen_q,1:sizesVector[chosen_q]])
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],feTemplatesMatrix[chosen_q,1:sizesVector[chosen_q]],lwd=1.5)
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],spectrumsMatrixEmissionLines[chosen_q,1:sizesVector[chosen_q]],col="gray45",lty=2)
}

#' @export
drawChosenPeaksComponents<-function(chosen_q, quasars, outputfit, wavelengthsMatrix, spectrumsMatrix, continuumsMatrix, feTemplatesMatrix, sizesVector){
  spectrumsMatrixORIG <- getSpectrumsMatrix(quasars)
  spectrumsMatrixNoIron <- rMinusMatrix(ptr, spectrumsMatrix, feTemplatesMatrix)
  spectrumsMatrixEmissionLines <- rMinusMatrix(ptr, spectrumsMatrixNoIron, continuumsMatrix)
  drawChosenSpectrumWithPeaksRawData(chosen_q, spectrumsMatrixEmissionLines, wavelengthsMatrix, outputfit, getParams(quasars), sizesVector)
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],spectrumsMatrix[chosen_q,1:sizesVector[chosen_q]],col="gray35")
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],continuumsMatrix[chosen_q,1:sizesVector[chosen_q]],col="gray45")
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],feTemplatesMatrix[chosen_q,1:sizesVector[chosen_q]],col="gray35")
}


#' @export
drawAllPeaksComponents2Files<-function(outpath, NSET, quasars, wavelengthsMatrix, spectrumsMatrix, continuumsMatrix, feTemplatesMatrix, sizesVector){
  spectrumsMatrixORIG <- getSpectrumsMatrix(quasars)
  spectrumsMatrixNoIron <- rMinusMatrix(ptr, spectrumsMatrix, feTemplatesMatrix)
  spectrumsMatrixEmissionLines <- rMinusMatrix(ptr, spectrumsMatrixNoIron, continuumsMatrix)
  qParams <- getParams(quasars)
  for(chosen_q in NSET){
    print("Zapisuję do jpg " %&% paste0("spSpec-",qParams[[chosen_q]]$mjd,"-",formatC(qParams[[chosen_q]]$plate, width=4, flag="0"),"-",formatC(qParams[[chosen_q]]$fiber, width=3, flag="0"),".fit")
    )
    jpeg(file=paste0(outpath,"quasar-",qParams[[chosen_q]]$mjd,"-",formatC(qParams[[chosen_q]]$plate, width=4, flag="0"),"-",formatC(qParams[[chosen_q]]$fiber, width=3, flag="0"),".jpg"), width = 500, height = 400, quality = 55, bg = "white")
  drawChosenSpectrumWithPeaksRawData(chosen_q, spectrumsMatrixEmissionLines, wavelengthsMatrix, outputfit, qParams, sizesVector)
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],spectrumsMatrixORIG[chosen_q,1:sizesVector[chosen_q]],lwd=.1,col="gray35")
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],spectrumsMatrix[chosen_q,1:sizesVector[chosen_q]])
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],continuumsMatrix[chosen_q,1:sizesVector[chosen_q]],col="gray45")
  lines(wavelengthsMatrix[chosen_q,1:sizesVector[chosen_q]],feTemplatesMatrix[chosen_q,1:sizesVector[chosen_q]],col="gray35")
    dev.off();
  }
}

#' @export
rContinuumTest <- function(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes,
						 windows, ampWavelength) 
{
	filteredMatrices <- rFilterWithWavelengthWindows(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes, windows)
		
	filteredMatrices <- rFilterNonpositive(quasarcl, filteredMatrices$spectrumsMatrix, 
										   filteredMatrices$wavelengthsMatrix, filteredMatrices$errorsMatrix, sizes)

	newSizes <- rCountIfNotInf(quasarcl, filteredMatrices$spectrumsMatrix)
	maxSize <- max(newSizes)
	
	spectrumsMatrixFiltered <- rCopyIfNotInf(quasarcl, filteredMatrices$spectrumsMatrix, maxSize)
	wavelengthsMatrixFiltered <- rCopyIfNotInf(quasarcl, filteredMatrices$aMatrix, maxSize)
	errorsMatrixFiltered <- rCopyIfNotInf(quasarcl, filteredMatrices$bMatrix, maxSize)
	
	spectrumsMatrixFilteredCopy <- rLog10(quasarcl, spectrumsMatrixFiltered)
	wavelengthsMatrixFilteredCopy <- rLog10(quasarcl, wavelengthsMatrixFiltered)
	
	cReglinResults <- rReglin(quasarcl, wavelengthsMatrixFilteredCopy, spectrumsMatrixFilteredCopy, newSizes)

	if(ampWavelength > (.Machine$double.xmin ^ 10.001))
	{
		lampLog10 <- log10(ampWavelength)
		wavelengthsMatrixFilteredCopy <- rMinusScalar(quasarcl, wavelengthsMatrixFilteredCopy, lampLog10)
	}

	reglinResults <- rReglin(quasarcl, wavelengthsMatrixFilteredCopy, spectrumsMatrixFilteredCopy, newSizes)
	reglin <- rFixReglinResults(quasarcl, cReglinResults, reglinResults)
	
	continuumMatrixFiltered<- rCalcCw(quasarcl, wavelengthsMatrixFiltered, reglin$cReglinResults)
	chisqsFiltered <- rChisq(quasarcl, spectrumsMatrixFiltered, continuumMatrixFiltered, errorsMatrixFiltered,
							 newSizes)
	chisqsFiltered <- rReduceContinuumChisqs(quasarcl, chisqsFiltered, newSizes)

	cfunDcfun <- rCalcCfunDcfun(quasarcl, wavelengthsMatrix, reglin$cReglinResults, reglin$reglinResults)

	return(list(dcontinuumsMatrix = cfunDcfun$dcontinuums, 
				continuumsMatrix = cfunDcfun$continuums,
				continuumChisqs = chisqsFiltered,
				continuumReglin = reglin$reglinResults,
				reglin = reglin$cReglinResults))
}

#' @export
rFeFitTest <- function(quasarcl, spectrumsData, feTemplate, feWindows, fitParameters)
{
	spectrumsMatrix <- spectrumsData$spectrumsMatrix
	wavelengthsMatrix <- spectrumsData$wavelengthsMatrix
	errorsMatrix <- spectrumsData$errorsMatrix
	continuumsMatrix <- spectrumsData$continuumsMatrix
	sizes <- spectrumsData$sizes
	
	templateFeMatrix <- rCalcFeTemplateMatrix(quasarcl, wavelengthsMatrix, sizes, feTemplate, fitParameters)
	templateFeMatrixCopy <- templateFeMatrix

	if (!(fitParameters$isSubC))
	{
		spectrumsMatrix <- rMinusMatrix(quasarcl, spectrumsMatrix, continuumsMatrix)
	}

	if (fitParameters$fitType == "WIN" || fitParameters$fitType == "FWIN")
	{
		filteredMatrices <- rFilterWithWavelengthWindows(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix,
														 sizes, feWindows) 
		spectrumsMatrix <- filteredMatrices$spectrumsMatrix
		wavelengthsMatrix <- filteredMatrices$wavelengthsMatrix
		errorsMatrix <- filteredMatrices$errorsMatrix
		
		filteredMatrices <- rFilterInfs(quasarcl, spectrumsMatrix, continuumsMatrix, templateFeMatrix, sizes)
		spectrumsMatrix <- filteredMatrices$spectrumsMatrix
		continuumsMatrix <- filteredMatrices$aMatrix
		templateFeMatrix <- filteredMatrices$bMatrix
	}

	if (fitParameters$fitType == "FWIN")
	{
		filteredMatrices <- rFilterZeros(quasarcl, templateFeMatrix, wavelengthsMatrix, errorsMatrix, sizes)
		templateFeMatrix <- filteredMatrices$spectrumsMatrix
		wavelengthsMatrix <- filteredMatrices$aMatrix
		errorsMatrix <- filteredMatrices$bMatrix
		
		filteredMatrices <- rFilterInf(quasarcl, templateFeMatrix, continuumsMatrix, spectrumsMatrix, sizes)
		templateFeMatrix <- filteredMatrices$spectrumsMatrix
		continuumsMatrix <- filteredMatrices$aMatrix
		spectrumsMatrix <- filteredMatrices$bMatrix
	}
	
	spectrumsMatrixFiltered <- spectrumsMatrix;
	wavelengthsMatrixFiltered <- wavelengthsMatrix;
	errorsMatrixFiltered <- errorsMatrix;
	continuumsMatrixFiltered <- continuumsMatrix;
	templateFeMatrixFiltered <- templateFeMatrix;
	sizesFeWindows <- sizes;
	maxSpectrumSize <- ASTRO_OBJ_SPEC_SIZE
	
	if (fitParameters$fitType == "WIN" || fitParameters$fitType == "FWIN")
	{
		sizesFeWindows <- rCountIfNotInf(quasarcl, filteredMatrices$spectrumsMatrix)
		maxSize <- max(sizesFeWindows)
		
		spectrumsMatrixFiltered <- rCopyIfNotInf(quasarcl, spectrumsMatrix, maxSize)
		wavelengthsMatrixFiltered <- rCopyIfNotInf(quasarcl, wavelengthsMatrix, maxSize)
		errorsMatrixFiltered <- rCopyIfNotInf(quasarcl, errorsMatrix, maxSize)
		continuumsMatrixFiltered <- rCopyIfNotInf(quasarcl, continuumsMatrix, maxSize)
		templateFeMatrixFiltered <- rCopyIfNotInf(quasarcl, templateFeMatrix, maxSize)
	}

	scaleRates <- rCalcFeTemplateScaleRates(quasarcl, spectrumsMatrixFiltered, 
												templateFeMatrixFiltered, sizesFeWindows,
												fitParameters)
	
	templateFeMatrixFiltered <- rMultiplyCol(quasarcl, templateFeMatrixFiltered, scaleRates)
	templateFeMatrixCopy <- rMultiplyCol(quasarcl, templateFeMatrixCopy, scaleRates)
	
	reducedChisqsFiltered <- rCalcReducedChisqs(quasarcl, spectrumsMatrixFiltered, templateFeMatrixFiltered,
												errorsMatrixFiltered, sizesFeWindows)
												
												
	#pobranie pelnych danych dla dalszych obliczeń
	spectrumsMatrix <- spectrumsData$spectrumsMatrix
	wavelengthsMatrix <- spectrumsData$wavelengthsMatrix
	errorsMatrix <- spectrumsData$errorsMatrix
	continuumsMatrix <- spectrumsData$continuumsMatrix
	
	if (!(fitParameters$isSubC))
	{
		spectrumsMatrix <- rMinusMatrix(quasarcl, spectrumsMatrix, continuumsMatrix)
	}
	
	reducedChisqsFull <- rCalcReducedChisqs(quasarcl, spectrumsMatrix, templateFeMatrixCopy,
											errorsMatrix, sizes)
		
	temp <- rDivide(quasarcl, templateFeMatrixCopy, continuumsMatrix)
	ewsFull <- rTrapz(quasarcl, temp, wavelengthsMatrix, sizes)
	
	templateFeMatrix <- templateFeMatrixCopy
	
	if (fitParameters$fitType == "WIN" || fitParameters$fitType == "FWIN")
	{
		filteredMatrices <- rFilterWithWavelengthWindows(quasarcl, spectrumsMatrix, wavelengthsMatrix, 
														 errorsMatrix, sizes, list(fitParameters$feFitRange))

		spectrumsMatrix <- filteredMatrices$spectrumsMatrix
		wavelengthsMatrix <- filteredMatrices$wavelengthsMatrix
		errorsMatrix <- filteredMatrices$errorsMatrix 
		
		filteredMatrices <- rFilterInfs(quasarcl, spectrumsMatrix, continuumsMatrix, templateFeMatrix, 
										sizes)
										
		templateFeMatrix <- filteredMatrices$bMatrix
		continuumsMatrix <- filteredMatrices$aMatrix
		spectrumsMatrix <- filteredMatrices$spectrumsMatrix
		
		sizesFiltered <- rCountIfNotInf(quasarcl, spectrumsMatrix)
		maxSize <- max(sizesFiltered)
		
		spectrumsMatrix <- rCopyIfNotInf(quasarcl, spectrumsMatrix, maxSize)
		wavelengthsMatrix <- rCopyIfNotInf(quasarcl, wavelengthsMatrix, maxSize)
		errorsMatrix<- rCopyIfNotInf(quasarcl, errorsMatrix, maxSize)
		continuumsMatrix <- rCopyIfNotInf(quasarcl, continuumsMatrix, maxSize)
		templateFeMatrix <- rCopyIfNotInf(quasarcl, templateFeMatrix, maxSize)
	
		reducedChisqsFeRange <- rCalcReducedChisqs(quasarcl, spectrumsMatrix, templateFeMatrix,
												   errorsMatrix, sizesFiltered)
												   
		templateFeMatrix <- rDivide(quasarcl, templateFeMatrix, continuumsMatrix)
		ewsFeRange <- rTrapz(quasarcl, templateFeMatrix, wavelengthsMatrix, sizesFiltered)
		
	} else 
	{
		reducedChisqsFeRange <- reducedChisqsFull
		ewsFeRange <- ewsFull
	}
	
	return (list(feTemplateMatrix = templateFeMatrixCopy, 
				 feScaleRates = scaleRates, 
				 feWindowsSizes = sizesFeWindows, 
				 feWindowsReducedChisqs = reducedChisqsFiltered,
				 feFullReducedChisqs = reducedChisqsFull, 
				 feFullEWs = ewsFull, 
				 feRangeReducedChisqs = reducedChisqsFeRange, 
				 feRangeEWs = ewsFeRange))
}

#' @export
rFitElementTest <- function(quasarcl, spectrumsLinesMatrix, continuumsMatrix, wavelengthsMatrix, 
							errorsMatrix, sizesVector, element)
{
	spectrumsLinesMatrixCopy <- spectrumsLinesMatrix
	wavelengthsMatrixCopy <- wavelengthsMatrix
	sizesCopy <- sizesVector
	
	filteredMatrices <- rFilterWithWavelengthWindows(quasarcl, spectrumsLinesMatrix, wavelengthsMatrix, 
													 continuumsMatrix, sizesVector, list(element$range))
	
	sizes <- rCountIfNotInf(quasarcl, filteredMatrices$spectrumsMatrix)
	maxSize <- max(sizes)

	spectrumsLinesMatrix<- rCopyIfNotInf(ptr, filteredMatrices$spectrumsMatrix, maxSize)
	wavelengthsMatrix <- rCopyIfNotInf(ptr, filteredMatrices$wavelengthsMatrix, maxSize)
	continuumsMatrix <- rCopyIfNotInf(ptr, filteredMatrices$errorsMatrix, maxSize)
	
	fitGResults <- (rep(list(element$fitGuess), nrow(continuumsMatrix)))
	fitGResults <- rFitGaussian(quasarcl, spectrumsLinesMatrix, wavelengthsMatrix, sizes, fitGResults)
	
	gaussiansMatrix <- rCalcGaussian(quasarcl, wavelengthsMatrix, fitGResults, sizes)
	gaussiansMatrix <- rDivide(quasarcl, gaussiansMatrix, continuumsMatrix)
	ews <- rTrapz(quasarcl, gaussiansMatrix, wavelengthsMatrix, sizes)

	gaussianChisqs <- rCalcGaussianChisqs(quasarcl, wavelengthsMatrixCopy, spectrumsLinesMatrixCopy,
										 errorsMatrix, fitGResults, sizesCopy)
		
	gaussianFWHMs <- rCalcGaussianFWHM(quasarcl, fitGResults)
	
	return(list(fitParams = fitGResults,
				ews = ews,
				chisqs = gaussianChisqs,
				gaussianFWHMs = gaussianFWHMs))
	
}

#' @export
rParameterizationTest <- function(ptr, quasars, options) 
{
	#domyślne wartości - funkcje zbiorcze rContinuum, rFeFit oraz rParameterization ładują je we własnym zakresie,
	#nie trzeba tego robić ręcznie, jeżeli z nich korzystamy
	
	spectralLines <- options$spectralLines
	continuumWindows <- options$continuumWindows
	ampWavelength <- options$ampWavelength
	feWindows<- options$feWindows
	fitParameters <- options$fitParameters	
	feTemplate <- options$feTemplate

	errorsMatrix <- getErrorsMatrix(quasars)
	
	sizesVector <- getSizesVector(quasars)  

	spectrumsMatrix <- getSpectrumsMatrix(quasars)
	spectrumsMatrix <- rCenteredMAVG(ptr, spectrumsMatrix, sizesVector, 25)
	
	
	abz <- getAbzParams(quasars)
	wavelengthsMatrix <- rGenerateWavelengthsMatrix(ptr, abz, ASTRO_OBJ_SPEC_SIZE)
	wavelengthsMatrix <- rTranspose(ptr, wavelengthsMatrix)

	
	filteredMatrices <- rFilterZeros(ptr, errorsMatrix, spectrumsMatrix, wavelengthsMatrix, sizesVector)
	newSizes <- rCountIfNotInf(ptr, filteredMatrices$bMatrix)
	maxSize <- max(newSizes)

	spectrumsMatrix<- rCopyIfNotInf(ptr, filteredMatrices$aMatrix, ASTRO_OBJ_SPEC_SIZE)
	wavelengthsMatrix <- rCopyIfNotInf(ptr, filteredMatrices$bMatrix, ASTRO_OBJ_SPEC_SIZE)
	errorsMatrix <- rCopyIfNotInf(ptr, filteredMatrices$spectrumsMatrix, ASTRO_OBJ_SPEC_SIZE)

	spectrumsMatrixCopy <- spectrumsMatrix
	
	
	max_iter <- 3
	for (i in seq(1:max_iter)) 
	{
		continuumResults <- rContinuumTest(ptr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix,
											newSizes, continuumWindows, ampWavelength)
		#continuumResults <- rContinuum(ptr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix,
		#								newSizes, continuumWindows, ampWavelength)

		spectrumsMatrix <- spectrumsMatrixCopy
	
		spectrumsData <- list(spectrumsMatrix = spectrumsMatrix,
							  wavelengthsMatrix = wavelengthsMatrix,
							  errorsMatrix = errorsMatrix,
							  continuumsMatrix = continuumResults$continuumsMatrix,
							  sizes = newSizes)

		feResults <- rFeFitTest(ptr, spectrumsData, feTemplate, feWindows, fitParameters)
		#feResults <- rFeFit(ptr, spectrumsData, feTemplate, feWindows, fitParameters)

		if (i < max_iter)
		{
			spectrumsMatrix <- rMinusMatrix(ptr, spectrumsMatrix, feResults$feTemplateMatrix)
		}
	}

	spectrumsEmissionLines <- rMinusMatrix(ptr, spectrumsMatrix, feResults$feTemplateMatrix)
	spectrumsEmissionLines <- rMinusMatrix(ptr, spectrumsEmissionLines, continuumResults$continuumsMatrix)
	
	fitElementsResults <- lapply(spectralLines, rFitElementTest, quasarcl = ptr, 
																 spectrumsLinesMatrix = spectrumsEmissionLines, 
																 continuumsMatrix = continuumResults$continuumsMatrix, 
																 wavelengthsMatrix = wavelengthsMatrix, 
																 errorsMatrix = errorsMatrix, 
																 sizesVector = newSizes)
	
	return (list(continuumChisqs = continuumResults$continuumChisqs,
				 continuumReglin = continuumResults$continuumReglin,
				 reglin = continuumResults$reglin,
				 feScaleRates = feResults$feScaleRates,
				 feWindowsSizes = feResults$feWindowsSizes,
				 feWindowsReducedChisqs = feResults$feWindowsReducedChisqs,
				 feFullReducedChisqs = feResults$feFullReducedChisqs,
				 feFullEWs = feResults$feFullEWs,
				 feRangeReducedChisqs = feResults$feRangeReducedChisqs,
				 feRangeEWs = feResults$feRangeEWs,
				 elementsFits = fitElementsResults,
				 spectrumsEmissionLines = spectrumsEmissionLines,
				 fitElementsResults = fitElementsResults,
				 abz = abz,
				 newSizes = newSizes,
				 wavelengthsMatrix = wavelengthsMatrix,
				 spectrumsMatrix=spectrumsMatrix
				 ))

}

