#ifndef QUASARCL_PARAMETERIZATION_H_
#define QUASARCL_PARAMETERIZATION_H_


#include "QuasarCL.hpp"

#include "quasarclTools.hpp"
#include "quasarclSpectrums.hpp"
#include "quasarclFeFit.hpp"
#include "quasarclContinuumRun.hpp"
#include "quasarclGaussian.hpp"
#include "quasarclBasics.hpp"
#include "quasarclMAVG.hpp"

#ifndef __CL_ENABLE_EXCEPTIONS
#define __CL_ENABLE_EXCEPTIONS
#endif

namespace quasarcl {

	
	Rcpp::List fitElement( 
			Rcpp::XPtr<QuasarCL> quasarclPtr,
			cl::Buffer& spectrumsEmissionLines, cl::Buffer& spectrumsEmissionLinesCopy,
			cl::Buffer& continuumsMatrix, /*cl::Buffer& dcontinuumsMatrix, błędy continuum nie są używane */
			cl::Buffer& wavelengthsMatrix, Rcpp::NumericMatrix wavelengthsMatrixHost,
			cl::Buffer& errorsMatrix,
			cl::Buffer& sizes, Rcpp::IntegerVector sizesHost, size_t size, Rcpp::List element)
	{
		auto queue = quasarclPtr->getQueue();
		auto context = quasarclPtr->getContext();
		
		// ________________________________________________________________________ //
		// __filtrowanie danych dla zakresu elementu_______________________________ //
		// ________________________________________________________________________ //
		cl_double2 range = Rcpp::as<cl_double2>(element["range"]);
		std::vector<cl_double2> feFitRangeVec(1, range);
		auto feFitRangeBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, feFitRangeVec.size() * sizeof(cl_double2));
		cl::copy(queue, feFitRangeVec.begin(), feFitRangeVec.end(), feFitRangeBuffer);

		filterWithWavelengthWindows(*quasarclPtr, spectrumsEmissionLines, wavelengthsMatrix, continuumsMatrix, sizes, size,
				feFitRangeBuffer, feFitRangeVec.size());
		// Zliczanie elementów po przefiltrowaniu.
		countIfNotInf(*quasarclPtr, spectrumsEmissionLines, size, ASTRO_OBJ_SPEC_SIZE, sizes);

		// ________________________________________________________________________ //
		// __Skracanie danych po filtracji_________________________________________ //
		// ________________________________________________________________________ //

		// Pobranie rozmiarów i wybranie największego (max), następnie
		// obliczenie najbliżej wielokrotności 64 większej od maksimum.
		// Potrzebne, żeby adresy elementów początkowych były wielokrotnością 64;
		Rcpp::IntegerVector sizes_filtered_vec(size, 0);
		cl::copy(queue, sizes, sizes_filtered_vec.begin(), sizes_filtered_vec.end());
		auto max = *std::max_element(sizes_filtered_vec.begin(), sizes_filtered_vec.end());

		size_t max_spectrum_filtered_size = max > 0 ? max : 64;
		max_spectrum_filtered_size = calcGlobalSize(64, max_spectrum_filtered_size);

		copyIfNotInf(*quasarclPtr, spectrumsEmissionLines, size, ASTRO_OBJ_SPEC_SIZE, spectrumsEmissionLines,
				max_spectrum_filtered_size);
		copyIfNotInf(*quasarclPtr, wavelengthsMatrix, size, ASTRO_OBJ_SPEC_SIZE, wavelengthsMatrix, max_spectrum_filtered_size);
		copyIfNotInf(*quasarclPtr, continuumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, continuumsMatrix, max_spectrum_filtered_size);
	//	tools.copyIfNotInf(dcontinuumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, dcontinuumsMatrix, max_spectrum_filtered_size);

		std::vector<cl_double4> fitGResultsVec(size);
		for (auto& i : fitGResultsVec)
		{
			i = element["fitGuess"];
		}

		auto fitGResults = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double4));
		cl::copy(queue, fitGResultsVec.begin(), fitGResultsVec.end(), fitGResults);

		fitGaussian(*quasarclPtr, spectrumsEmissionLines, wavelengthsMatrix, sizes, size, fitGResults);
		cl::copy(queue, fitGResults, fitGResultsVec.begin(), fitGResultsVec.end());

		// spectrumsEmissionLines już nie jest potrzebne
		auto gaussiansMatrix = spectrumsEmissionLines;
		calcGaussian(*quasarclPtr, wavelengthsMatrix, fitGResults, sizes, size, max_spectrum_filtered_size, gaussiansMatrix);
		
		// ________________________________________________________________________ //
		// __Gorne oszacowanie błędu na EW korzystając różniczki zupelnej__________ //
		// ________________________________________________________________________ //

		Rcpp::NumericVector ewsVec(size);
		{
			auto ews = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double));;
			divide(*quasarclPtr, gaussiansMatrix, size, max_spectrum_filtered_size, continuumsMatrix, gaussiansMatrix);
			trapz(*quasarclPtr, gaussiansMatrix, wavelengthsMatrix, size, max_spectrum_filtered_size, sizes, ews);
			cl::copy(queue, ews, ewsVec.begin(), ewsVec.end());
		}

		// ________________________________________________________________________ //
		// __Dopasowanie do Gaussiana______________________________________________ //
		// ________________________________________________________________________ //

		// Kopiowanie pełnych długości fal i rozmiarów
		cl::copy(queue, wavelengthsMatrixHost.begin(), wavelengthsMatrixHost.end(), wavelengthsMatrix);
		cl::copy(queue, sizesHost.begin(), sizesHost.end(), sizes);

		// Obliczenie chi kwadrat
		Rcpp::NumericVector gaussianChisqsVec(size);
		{
			auto gaussianChisqs = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double));
			calcGaussianChisq(*quasarclPtr, wavelengthsMatrix, spectrumsEmissionLinesCopy, errorsMatrix, fitGResults, sizes, size,
					gaussianChisqs);
			cl::copy(queue, gaussianChisqs, gaussianChisqsVec.begin(), gaussianChisqsVec.end());
		}
		// ________________________________________________________________________ //
		// __Obliczenie FWHM_______________________________________________________ //
		// ________________________________________________________________________ //

		Rcpp::NumericVector gaussianFWHMsVec(size);
		{
			// Obliczamy
			auto gaussianFWHMs = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double));;
			calcGaussianFWHM(*quasarclPtr, fitGResults, gaussianFWHMs, size);
			cl::copy(queue, gaussianFWHMs, gaussianFWHMsVec.begin(), gaussianFWHMsVec.end());
		}
		
		return Rcpp::List::create(Rcpp::Named("fitParams") = fitGResultsVec,
					Rcpp::Named("ews") = ewsVec,
					Rcpp::Named("chisqs") = gaussianChisqsVec,
					Rcpp::Named("gaussianFWHMs") = gaussianFWHMsVec
		);
	}

	
	
	
	Rcpp::List parameterization(Rcpp::XPtr<QuasarCL> quasarclPtr, Rcpp::NumericMatrix spectrumsMatrixHost, 
								Rcpp::NumericMatrix errorsMatrixHost,
								Rcpp::IntegerVector sizesHost, std::vector<cl_double4> abz_, 
								Rcpp::List elements,
								std::vector<cl_double2> continuumWindows, cl_double ampWavelength,
								std::vector<cl_double2> feWindows, 
								Rcpp::List feTemplate, Rcpp::List fitParameteres)
	{
		auto queue = quasarclPtr->getQueue();
		auto context = quasarclPtr->getContext();
		
		size_t size = spectrumsMatrixHost.rows();
		// ________________________________________________________________________ //
		// __bufory ze danymi______________________________________________________ //
		// ________________________________________________________________________ //

		cl::Buffer spectrumsMatrix = cl::Buffer(context, CL_MEM_READ_WRITE, ASTRO_OBJ_SPEC_SIZE * size * sizeof(cl_double));
		cl::Buffer wavelengthsMatrix = cl::Buffer(context, CL_MEM_READ_WRITE, ASTRO_OBJ_SPEC_SIZE * size * sizeof(cl_double));
		cl::Buffer errorsMatrix = cl::Buffer(context, CL_MEM_READ_WRITE, ASTRO_OBJ_SPEC_SIZE * size * sizeof(cl_double));
		cl::Buffer sizes = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_uint));
		cl::Buffer continuumsMatrix;
		cl::Buffer feTemplatesMatrix;

		// ________________________________________________________________________ //
		// __bufory ze danymi______________________________________________________ //
		// ________________________________________________________________________ //

		Rcpp::NumericMatrix wavelengthsMatrixHost(size, ASTRO_OBJ_SPEC_SIZE);
		Rcpp::NumericMatrix continuumsMatrixHost(size, ASTRO_OBJ_SPEC_SIZE);
		Rcpp::NumericMatrix dcontinuumsMatrixHost(size, ASTRO_OBJ_SPEC_SIZE);
		Rcpp::NumericMatrix feTemplatesMatrixHost(size, ASTRO_OBJ_SPEC_SIZE);
	
		//
		// __generowanie długości fal______________________________________________ //
		// ________________________________________________________________________ //
		{
			auto abz_buffer =  cl::Buffer(context, CL_MEM_READ_ONLY, size * sizeof(cl_double4));
			auto temp =  cl::Buffer(context, CL_MEM_READ_WRITE, ASTRO_OBJ_SPEC_SIZE * size * sizeof(cl_double));

			cl::copy(queue, abz_.begin(), abz_.end(), abz_buffer);
			generateWavelengthsMatrix(*quasarclPtr, abz_buffer, size, temp);
			transpose(*quasarclPtr, temp, ASTRO_OBJ_SPEC_SIZE, size, wavelengthsMatrix);
			
		}
		
		// ________________________________________________________________________ //
		// __kopiowanie danych do buforów__________________________________________ //
		// ________________________________________________________________________ //

		cl::copy(queue, spectrumsMatrixHost.begin(), spectrumsMatrixHost.end(), spectrumsMatrix);
		cl::copy(queue, errorsMatrixHost.begin(), errorsMatrixHost.end(), errorsMatrix);
		cl::copy(queue, sizesHost.begin(), sizesHost.end(), sizes);
		
		
		// ________________________________________________________________________ //
		// __transpozycja danych (i wygładzenie)___________________________________ //
		//___nie ma potrzeby transponować danych, są we właściwym ułożeniu_________ //
		// ________________________________________________________________________ //
		{
			auto temp = cl::Buffer(context, CL_MEM_READ_WRITE, ASTRO_OBJ_SPEC_SIZE * size * sizeof(cl_double));		
			centeredMAVG(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, sizes, temp, 10);
			spectrumsMatrix = temp;
			
			
		}
		// ________________________________________________________________________ //
		// __filtrowanie wartości, dla których błąd wynosi 0_______________________ //
		// ________________________________________________________________________ //
		{
			// filtrowanie
			
			
			filterZeros(*quasarclPtr, errorsMatrix, spectrumsMatrix, wavelengthsMatrix, sizes, size);
			// ustalenie nowych rozmiarów
			countIfNotInf(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, sizes);
			
			// kopiowanie przefiltrowanych danych
			copyIfNotInf(*quasarclPtr, errorsMatrix, size, ASTRO_OBJ_SPEC_SIZE, errorsMatrix, ASTRO_OBJ_SPEC_SIZE);	
			copyIfNotInf(*quasarclPtr, wavelengthsMatrix, size, ASTRO_OBJ_SPEC_SIZE, wavelengthsMatrix, ASTRO_OBJ_SPEC_SIZE);
			copyIfNotInf(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, spectrumsMatrix, ASTRO_OBJ_SPEC_SIZE);
		}
		
		// ________________________________________________________________________ //
		// __kopiowanie transponowanych, wygładzonych, przefiltrowanych danych_____ //
		// __na pamięć hosta_______________________________________________________ //
		// ________________________________________________________________________ //
		
		cl::copy(queue, sizes, sizesHost.begin(), sizesHost.end());
		cl::copy(queue, spectrumsMatrix, spectrumsMatrixHost.begin(), spectrumsMatrixHost.end());
		cl::copy(queue, wavelengthsMatrix, wavelengthsMatrixHost.begin(), wavelengthsMatrixHost.end());
		cl::copy(queue, errorsMatrix, errorsMatrixHost.begin(), errorsMatrixHost.end());
		
		
		ContinuumResults continuumResults;
		FitResults feResults;

		//
		int max_i = 3;
		for (int i = 0; i < max_i; i++)
		{

			// ________________________________________________________________________ //
			// __obliczenie kontinuum__________________________________________________ //
			// ________________________________________________________________________ //
			continuumResults = continuum(quasarclPtr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes, size, ASTRO_OBJ_SPEC_SIZE, continuumWindows,
					ampWavelength);

			continuumsMatrix = continuumResults.continuumsMatrix;
			
			
			// ________________________________________________________________________ //
			// __kopiowanie kontinuum na host__________________________________________ //
			// ________________________________________________________________________ //

			cl::copy(queue, continuumsMatrix, continuumsMatrixHost.begin(), continuumsMatrixHost.end());
			
			// Błędy kontinuum NIE są UŻYWANE
			//
			// Jeżeli ostatnia pętla to trzeba zapisać błędy kontinuum na host
	//		if (i == (max_i - 1))
	//		{
	//			cl::copy(queue, continuumResult.dcontinuumsMatrix,
	//					dcontinuumsMatrixHost.begin(), dcontinuumsMatrixHost.end());
	//		}

			// ________________________________________________________________________ //
			// __wczytywanie do buforów oryginalnych danych____________________________ //
			// ________________________________________________________________________ //

			// Trzeba zadeklarować nowy bufor, bo poprzedni zajęła macierz continuumsMatrix
			spectrumsMatrix = cl::Buffer(context, CL_MEM_READ_WRITE, size * ASTRO_OBJ_SPEC_SIZE * sizeof(cl_double));

			cl::copy(queue, spectrumsMatrixHost.begin(), spectrumsMatrixHost.end(), spectrumsMatrix);
			cl::copy(queue, wavelengthsMatrixHost.begin(), wavelengthsMatrixHost.end(), wavelengthsMatrix);
			cl::copy(queue, errorsMatrixHost.begin(), errorsMatrixHost.end(), errorsMatrix);
			cl::copy(queue, sizesHost.begin(), sizesHost.end(), sizes);

			// ________________________________________________________________________ //
			// __dopasowanie do szablonu żelaza________________________________________ //
			// ________________________________________________________________________ //


			/*Data specData = {
				spectrumsMatrixHost, wavelengthsMatrixHost,
				errorsMatrixHost, continuumsMatrixHost, sizesHost
			};*/

			Rcpp::List specData = Rcpp::List::create(Named("spectrumsMatrix") = spectrumsMatrixHost,
													 Named("wavelengthsMatrix") = wavelengthsMatrixHost,
													 Named("errorsMatrix") = errorsMatrixHost,
													 Named("continuumsMatrix") = continuumsMatrixHost,
													 Named("sizes") = sizesHost);
			
			Buffers specBuffers = {
					spectrumsMatrix, wavelengthsMatrix,
					errorsMatrix, continuumsMatrix, sizes
			};

			feResults = feFit(quasarclPtr, specData, specBuffers, size, feTemplate, feWindows, fitParameteres);

			if (i < (max_i - 1))
			{
				cl::copy(queue, spectrumsMatrixHost.begin(), spectrumsMatrixHost.end(), spectrumsMatrix);
				cl::copy(queue, wavelengthsMatrixHost.begin(), wavelengthsMatrixHost.end(), wavelengthsMatrix);
				cl::copy(queue, errorsMatrixHost.begin(), errorsMatrixHost.end(), errorsMatrix);
				cl::copy(queue, sizesHost.begin(), sizesHost.end(), sizes);

				minus(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, feResults.feTemplateMatrix, spectrumsMatrix);
			}

			// Jeżeli ostatnia pętla to trzeba zapisać feTemplateMatrix
			if (i == (max_i - 1))
			{
				cl::copy(queue, feResults.feTemplateMatrix, feTemplatesMatrixHost.begin(), feTemplatesMatrixHost.end());
			}
		}

		// Błędy kontinuum NIE są UŻYWANE
	//	auto dcontinuumsMatrix = utils::make_qbuffer(context, CL_MEM_READ_WRITE, size);
		cl::copy(queue, continuumsMatrixHost.begin(), continuumsMatrixHost.end(), continuumsMatrix);
		// Błędy kontinuum NIE są UŻYWANE
	//	cl::copy(queue, dcontinuumsMatrixHost.begin(), dcontinuumsMatrixHost.end(), dcontinuumsMatrix);
		cl::copy(queue, spectrumsMatrixHost.begin(), spectrumsMatrixHost.end(), spectrumsMatrix);
		cl::copy(queue, wavelengthsMatrixHost.begin(), wavelengthsMatrixHost.end(), wavelengthsMatrix);
		cl::copy(queue, errorsMatrixHost.begin(), errorsMatrixHost.end(), errorsMatrix);
		cl::copy(queue, sizesHost.begin(), sizesHost.end(), sizes);

		// ________________________________________________________________________ //
		// __obliczenie linii emisyjnych___________________________________________ //
		// ________________________________________________________________________ //

		auto spectrumsEmissionLines = feResults.feTemplateMatrix;
		
		minus(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, feResults.feTemplateMatrix, spectrumsEmissionLines);
		minus(*quasarclPtr, spectrumsEmissionLines, size, ASTRO_OBJ_SPEC_SIZE, continuumsMatrix, spectrumsEmissionLines);

		// ________________________________________________________________________ //
		// __zapisanie uzyskanych linii emisyjnych_________________________________ //
		// ________________________________________________________________________ //

		auto spectrumsEmissionLinesCopy = spectrumsMatrix;
		queue.enqueueCopyBuffer(spectrumsEmissionLines, spectrumsEmissionLinesCopy, 0, 0,
				sizeof(cl_double) * size * ASTRO_OBJ_SPEC_SIZE);

		// Dopasowanie do elementów
		Rcpp::List fitElementsResults;
		for (auto& element : elements)
		{
			auto fit = fitElement(quasarclPtr, spectrumsEmissionLines, spectrumsEmissionLinesCopy,
					continuumsMatrix, /* dcontinuumsMatrix, */
					wavelengthsMatrix, wavelengthsMatrixHost, errorsMatrix,
					sizes, sizesHost, size, element);
			fitElementsResults.push_back(fit);
			
			// Kopiowanie pełnych danych, żeby można ich użyć dla kolejnego elementu
			// w kolejnej iteracji pętli.

			queue.enqueueCopyBuffer(spectrumsEmissionLinesCopy, spectrumsEmissionLines, 0, 0,
					sizeof(cl_double) * size * ASTRO_OBJ_SPEC_SIZE);
			cl::copy(queue, continuumsMatrixHost.begin(), continuumsMatrixHost.end(), continuumsMatrix);

		}
		
		Rcpp::List results = Rcpp::List::create(Rcpp::Named("continuumChisqs") = continuumResults.chisqs,
				Rcpp::Named("continuumReglin") = continuumResults.c_reglin_result,
				Rcpp::Named("reglin") = continuumResults.reglin_result,
				Rcpp::Named("feScaleRates") = feResults.scaleRates,
				Rcpp::Named("feWindowsSizes") = feResults.sizes_fewindows,
				Rcpp::Named("feWindowsReducedChisqs") = feResults.reducedChisqs_fewindows,
				Rcpp::Named("feFullReducedChisqs") = feResults.reducedChisqs_full,
				Rcpp::Named("feFullEWs") = feResults.ews_full,
				Rcpp::Named("feRangeReducedChisqs") = feResults.reducedChisqs_feRange,
				Rcpp::Named("feRangeEWs") = feResults.ews_feRange,
				Rcpp::Named("elementsFits") = fitElementsResults,
				Rcpp::Named("wavelengthsMatrix") = wavelengthsMatrixHost,
				Rcpp::Named("spectrumsMatrix") = spectrumsMatrixHost,
				Rcpp::Named("continuumsMatrix") = continuumsMatrixHost,
				Rcpp::Named("feTemplatesMatrix") = feTemplatesMatrixHost,
				Rcpp::Named("sizes") = sizesHost);
		return results;
	}
	
}

#endif
