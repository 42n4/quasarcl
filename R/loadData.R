#' @importFrom utils file_test


#' @export
loadSpectralLines <- function(file)
{
	if(!file_test("-f", file))
	{
        stop("File with spectral lines does not exist")
    }
	df <- read.table(file, header = FALSE)
	colnames(df) <- c("name", "range0", "range1", "fitGuess0", "fitGuess1", "fitGuess2" )
	spectralLines <- apply(df[, 1:length(df)], 1, specLineToList)
	return (spectralLines)
}



#' @export
loadIronEmissionTemplate <- function(file) 
{
	if(!file_test("-f", file))
	{
        stop("File with template does not exist")
    }
	template <- read.table(file, header = FALSE)
	colnames(template) <- c("wavelengths", "values")

	return (list(wavelengths = template[, "wavelengths"], 
				 values = template[, "values"]))
}



#' @export
loadWindows <-function(file) 
{
	if(!file_test("-f", file))
	{
        stop("File with  windows does not exist")
    }
    df <- read.table(file, header = FALSE)
    return (split(as.matrix(df), seq(nrow(df))))
}

#' @export
loadQuasarsFromFitFiles <- function(path, from = 1, N = 10000) 
{
	fitFiles <- list.files(path=path, pattern="\\.fit$", full.names = TRUE)
	if(length(fitFiles) < from) stop("za mało plików")
	if(length(fitFiles) < (from+N-1)) N <- length(fitFiles)-from+1
	fitFiles <- fitFiles[from:(from+N-1)]
	print(paste0("Wczytuję ", N, " kwazarów z fit"));
	if (identical(fitFiles, character(0))) 
	{
		stop("No .fit files found")
	}
	return (lapply(fitFiles, parseFitFile))
}

#' @export
loadQuasarsFromSelectedFitFiles <- function(fitFiles, from = 1, N = 10000) 
{
	if(length(fitFiles) < from) stop("za mało plików")
	if(length(fitFiles) < (from+N-1)) N <- length(fitFiles)-from+1
	fitFiles <- fitFiles[from:(from+N-1)]
	print(paste0("Wczytuję ", N, " kwazarów z fit"));
	if (identical(fitFiles, character(0))) 
	{
		stop("No .fit files found")
	}
	return (lapply(fitFiles, parseFitFile))
}

#' @export
downloadFitFilesFromServer <- function(path, from = 1, N =-1)
{
	if (!dir.exists(path)) 
	{
		stop("Directory does not exist")
	}
	df <-read.table(SERVER_QUASARS_INFO, nrows=N, fill=TRUE , na.strings="", skip=(80+from-1))
	apply(df, 1, downloadFitFile, path)
}



#' @export
getSpectrumsMatrix <- function(quasars) 
{
	spectrums <- lapply(lapply(quasars,  `[[`, 'values' ), 
						function(spectrum) { c(spectrum, rep(Inf, ASTRO_OBJ_SIZE - length(spectrum))) })
	return(do.call(rbind, spectrums))
}



#' @export
getErrorsMatrix <- function(quasars) 
{
	errors <- lapply(lapply(quasars,  `[[`, 'error' ), 
					 function(error) { c(error, rep(Inf, ASTRO_OBJ_SIZE - length(error))) })
	return(do.call(rbind, errors))
}



#' @export
getSizesVector <-function(quasars) 
{
	sizes <-sapply(lapply(quasars,  `[[`, 'values' ), length)
	return (sizes)
}



#' @export
getAbzParams <- function(quasars) 
{
	abz <- lapply(getParams(quasars), function(param) { c(param$a, param$b, param$z, 0)})
	return(abz)
}


#' @export
getParams <- function(quasars) 
{
	params <-lapply(quasars, `[[`, 'params')
	return(params)
}



#' @export
getLambda <- function(path, z) 
{
	data<-"lambda.csv"
	assign(data,read.csv(paste("file://", path,data,sep=""),head=FALSE,sep=",",dec=".",na.strings=c("NA", "BD", "bd", "", "?")))
	lambda<-10^get(data)[,1]
	return(lambda/(z+1))
}
