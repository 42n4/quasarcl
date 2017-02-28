#include "RcppExtended.hpp"
#include "quasarclFeWindows.hpp"
#include "quasarclFeFit.hpp"
		
	

	
//[[Rcpp::export]]
SEXP cppFeFit(SEXP quasarclPtr_, SEXP spectrumsMatrix_, SEXP wavelengthsMatrix_, 
				 SEXP errorsMatrix_, SEXP continuumsMatrix_, SEXP sizesVector_, 
				 SEXP feTemplateList_, SEXP windowsVector_, SEXP fitParametersList_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
	Rcpp::NumericMatrix wavelengthsMatrix(wavelengthsMatrix_);
	Rcpp::NumericMatrix errorsMatrix(errorsMatrix_);
	Rcpp::NumericMatrix continuumsMatrix(errorsMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	Rcpp::List feTemplate(feTemplateList_);
	std::vector<cl_double2> windowsVector = Rcpp::as<std::vector<cl_double2> >(windowsVector_);
	Rcpp::List fitParameters(fitParametersList_);
	
	
	auto spectrumsNumber = spectrumsMatrix.rows();
	auto spectrumsSize = spectrumsMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferSpectrums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferWavelengths = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferErrors = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferContinuums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	
	cl::copy(queue, spectrumsMatrix.begin(), spectrumsMatrix.end(), bufferSpectrums);
	cl::copy(queue, wavelengthsMatrix.begin(), wavelengthsMatrix.end(), bufferWavelengths);
	cl::copy(queue, errorsMatrix.begin(), errorsMatrix.end(), bufferErrors);
	cl::copy(queue, continuumsMatrix.begin(), continuumsMatrix.end(), bufferContinuums);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	
	quasarcl::Data specData = {
		spectrumsMatrix,
		wavelengthsMatrix,
		errorsMatrix,
		continuumsMatrix,
		sizesVector
	};

	quasarcl::Buffers specBuffers = {
		bufferSpectrums,
		bufferWavelengths,
		bufferErrors,
		bufferContinuums,
		bufferSizes
	};

	auto results = quasarcl::feFit(quasarclPtr, specData, specBuffers, spectrumsNumber, feTemplate, 
									  windowsVector, fitParameters); 
	
	Rcpp::NumericMatrix templateFeMatrix(spectrumsNumber, spectrumsSize);
	cl::copy(queue, results.feTemplateMatrix, templateFeMatrix.begin(), templateFeMatrix.end());
	
	return Rcpp::List::create(Rcpp::Named("feTemplateMatrix") = templateFeMatrix,
							  Rcpp::Named("scaleRates") = results.scaleRates,
							  Rcpp::Named("sizes_fewindows") = results.sizes_fewindows,
							  Rcpp::Named("reducedChisqs_fewindows") = results.reducedChisqs_fewindows,
							  Rcpp::Named("reducedChisqs_full") = results.reducedChisqs_full,
							  Rcpp::Named("ews_full") = results.ews_full,
							  Rcpp::Named("reducedChisqs_feRange") = results.reducedChisqs_feRange,
							  Rcpp::Named("ews_feRange") = results.ews_feRange);
}



//[[Rcpp::export]]
SEXP cppCalcFeTemplateMatrix(SEXP quasarclPtr_, SEXP wavelengthsMatrix_, SEXP sizesVector_,
							 SEXP feTemplateList_, SEXP fitParametersList_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix wavelengthsMatrix(wavelengthsMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	Rcpp::List feTemplate(feTemplateList_);
	Rcpp::List fitParameters(fitParametersList_);

	
	auto spectrumsNumber = wavelengthsMatrix.rows();
	auto spectrumsSize = wavelengthsMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferWavelengths = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	
	cl::copy(queue, wavelengthsMatrix.begin(), wavelengthsMatrix.end(), bufferWavelengths);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	auto bufferOutput = quasarcl::calcFeTemplateMatrix(quasarclPtr, bufferWavelengths, bufferSizes, 
													   spectrumsNumber, feTemplate, fitParameters); 
	
	Rcpp::NumericMatrix outputMatrix(spectrumsNumber, spectrumsSize);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}



//[[Rcpp::export]]
SEXP cppCalcFeTemplateScaleRates(SEXP quasarclPtr_, SEXP spectrumsMatrix_, SEXP templateFeMatrix_,
								 SEXP sizesVector_, SEXP fitParametersList_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
	Rcpp::NumericMatrix templateFeMatrix(templateFeMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	Rcpp::List fitParameters(fitParametersList_);
	
	
	auto spectrumsNumber = spectrumsMatrix.rows();
	auto spectrumsSize = spectrumsMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferSpectrums = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferTemplateFe = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	
	cl::copy(queue, spectrumsMatrix.begin(), spectrumsMatrix.end(), bufferSpectrums);
	cl::copy(queue, templateFeMatrix.begin(), templateFeMatrix.end(), bufferTemplateFe);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	auto bufferOutput = quasarcl::calcFeTemplateScaleRates(quasarclPtr, bufferSpectrums, 
														   bufferTemplateFe, bufferSizes, spectrumsNumber, 
														   spectrumsSize, fitParameters); 
	
	Rcpp::NumericVector outputVector(spectrumsNumber);
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());	
	return outputVector;
}



//[[Rcpp::export]]
SEXP cppCalcReducedChisqs(SEXP quasarclPtr_, SEXP fMatrix_, SEXP yMatrix_, SEXP errorsMatrix_,
						  SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix fMatrix(fMatrix_);
	Rcpp::NumericMatrix yMatrix(yMatrix_);
	Rcpp::NumericMatrix errorsMatrix(errorsMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = fMatrix.rows();
	auto spectrumsSize = fMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferF = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferY = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferErrors = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	
	cl::copy(queue, fMatrix.begin(), fMatrix.end(), bufferF);
	cl::copy(queue, yMatrix.begin(), yMatrix.end(), bufferY);
	cl::copy(queue, errorsMatrix.begin(), errorsMatrix.end(), bufferErrors);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	auto bufferOutput = quasarcl::calcReducedChisqs(quasarclPtr, bufferF, bufferY, bufferErrors, 
													bufferSizes, spectrumsNumber, spectrumsSize); 
	
	Rcpp::NumericVector outputVector(spectrumsNumber);
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return outputVector;
}



//[[Rcpp::export]]
SEXP cppCpuConvolve(SEXP quasarclPtr_, SEXP signalVector_, SEXP kernelVector_, SEXP same_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericVector signalVector(signalVector_);
	Rcpp::NumericVector kernelVector(kernelVector_);
	bool same = Rcpp::as<bool>(same_);
	
	return quasarcl::cpuConvolve(signalVector, kernelVector, same); 
}



//[[Rcpp::export]]
SEXP cppReduceFeChisq(SEXP quasarclPtr_, SEXP chisqsVector_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericVector chisqsVector(chisqsVector_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto size = chisqsVector.size();
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferChisqs = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, size * sizeof(uint));
	
	cl::copy(queue, chisqsVector.begin(), chisqsVector.end(), bufferChisqs);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	quasarcl::reduceFeChisqs(*quasarclPtr, bufferChisqs, bufferSizes, size);
	
	Rcpp::NumericVector outputChisqsVector(size);
	cl::copy(queue, bufferChisqs, outputChisqsVector.begin(), outputChisqsVector.end());
	return outputChisqsVector;
}  