// This file is part of the reference implementation for the paper
//    Fast Lossy Compression of 3D Unit Vector Sets
//    Sylvain Rousseau and Tamy Boubekeur
//    SIGGRAPH Asia 2017 Technical Briefs (SA '17)
//    DOI: https://doi.org/10.1145/3145749.3149436
//
// All rights reserved. Use of this source code is governed by a 
// MIT license that can be found in the LICENSE file.


#include <iostream>
#include <ctime>
#include <vector>


#include "uniquant.h"
#include "sfibonacci.h"
#include "timer.h"

using CompressedVector = uint16_t;
using CompressedWindow = std::vector<CompressedVector>;

// since the inverse mapping have some numerical precision issue (see the original paper: Spherical Fibonacci Mapping,Keinert & al 2015) , 
// it's advisable to keep the number of fibonacci points under 2^22
const uint32_t nbFiboPoints = std::numeric_limits<CompressedVector>::max(); // should be available on the master and the workers
const size_t nbVectors = 1'000'000;
const unsigned level = 13; // should be available on the master and the workers



/// The source code in the utl namespace is not mandatory for the compression
namespace utl
{
    /// generate a random uniform unit vector using rejection method
	glm::fvec3 randomUnitVector()
	{
		glm::fvec3 randVec;
		do
		{
			randVec = 2.f * glm::fvec3(float(rand()) / RAND_MAX,
									   float(rand()) / RAND_MAX,
									   float(rand()) / RAND_MAX)
				- glm::fvec3(1.f, 1.f, 1.f);
		} while (glm::dot(randVec, randVec) > 1.f);
		return glm::normalize(randVec);
	}
	
	/// generate $nbVectors uniform random vectors using rejection method
	std::vector<glm::fvec3> randomUnitVectorsArray(const size_t nbVectors)
	{
		std::vector<glm::fvec3> data;
		data.reserve(nbVectors);
		for (unsigned i = 0; i < nbVectors; ++i)
			data.push_back(utl::randomUnitVector());
		return data;
	}

	/// unsigned angle in radians between $v1 and $v2
	double computeError(const glm::fvec3 & v1, const glm::fvec3 & v2)
	{
		glm::dvec3 v1d(v1), v2d(v2);
		v1d = glm::normalize(v1d);
		v2d = glm::normalize(v2d);
		double dp = glm::dot(v1d, v2d);
		// numerical imprecision, avoid nans
		if (dp > 1.0)
		dp = 1.0;
		return std::abs(std::acos(dp));
	}

	/// convert $radians from radian to degree
	double radians2degrees(const double & radians)
	{
		return radians * 180.0 / uvc::M_PI;
	}

}

/// Sexample of code that can be executed on a worker node, decompress unit vectors and return them
std::vector<std::vector<glm::fvec3>> sendToWorker(const std::vector<CompressedWindow> & compressedData)
{
	const uint32_t mask = uvc::NMostSignificantBitsMask(level);

	std::vector<std::vector<glm::fvec3>> uncompressedData;
	uncompressedData.resize(compressedData.size());
	for(int w = 0; w < int(compressedData.size()); ++w)
	{
		///compute average from w
		glm::fvec3 average = uvc::computeAverage(w, mask, level);
		glm::dvec3 daverage = glm::normalize(glm::dvec3(average));

		///uncompress unit vectors
		uncompressedData[w].resize(compressedData[w].size());
		

		const float ratio = uvc::computeRatio(average, uvc::NMostSignificantBitsMask(level));

		#pragma omp parallel for
		for (int i = 0; i < int(compressedData[w].size()); ++i)
		{
			glm::fvec3 v = fib::SF(compressedData[w][i], nbFiboPoints);
			if (abs(glm::dot(glm::normalize(glm::dvec3(v)), daverage)) < 1.0) // numerical instabilities
				v = uvc::uniformMapping(v, average, ratio);
			uncompressedData[w][i] = v;
		}
	}

	/// here you can work with your vectors, here, we are just sending back the uncompressed vectors to compute error,
	/// in a distributed rendering engine you could for instance send back the result of the intersection 
	return uncompressedData;
}

int main()
{
	
	Timer timer;
	srand(unsigned(time(NULL)));

	///////////////////////
	/// initialize data////
	///////////////////////
	std::cout << "Creating data \t\t\t...";
	timer.reset();
	std::vector<glm::fvec3> data = utl::randomUnitVectorsArray(nbVectors);
	timer.printElapsed(" done");

	////////////////////////
	/// group unit vectors//
	////////////////////////
	std::cout << "Grouping data \t\t\t...";
	timer.reset();
	std::vector<std::vector<size_t>> indices;
	std::vector<std::vector<glm::fvec3>> groupedData;

	uvc::group(data, indices, groupedData, level);
	timer.printElapsed(" done");

	/////////////////////
	/// compress them ///
	/////////////////////
	std::cout << "Compressing data \t\t...";
	timer.reset();
	std::vector<CompressedWindow> compressedData;
	compressedData.resize(groupedData.size());

	uint32_t mask = uvc::NMostSignificantBitsMask(level);


    #pragma omp parallel for
	for (int i = 0; i < int(groupedData.size()); ++i) // each window is processed in parallel
	{
		/// empty windows
		if (groupedData[i].size() == 0)
			continue;


		glm::fvec3 average = uvc::computeAverage(groupedData[i][0], mask);
		float ratio = uvc::computeRatio(average, mask);
		
		CompressedWindow compressedWin;
		compressedWin.resize(groupedData[i].size());

		glm::dvec3 daverage = glm::normalize(glm::dvec3(average));

		for (int j = 0; j < int(groupedData[i].size()); ++j)
		{
			if (abs(glm::dot(glm::normalize(glm::dvec3(groupedData[i][j])), daverage)) < 1.0) /// numerical issue
				groupedData[i][j] = glm::normalize(uvc::inverseUniformMapping(groupedData[i][j], average, ratio)); /// mapping
			compressedWin[j] = CompressedVector(fib::inverseSF(groupedData[i][j], nbFiboPoints)); ///quantization
		}
		compressedData[i] = compressedWin;
	}

	timer.printElapsed(" done");
	///////////////////////////////////////////////////////////////////////////////////////////
	//// END OF COMPRESSION ///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////
	/// send data to the worker //
	//////////////////////////////
	std::cout << "Decompressing data (workers) \t...";
	timer.reset();
	std::vector<std::vector<glm::fvec3>> result = sendToWorker(compressedData);
	timer.printElapsed(" done");
	


	///-------------------------------------------------------------------------
	/// The result of the calculation are returned by the worker, using indices  
	/// you can register the result in the correct place.
	///-------------------------------------------------------------------------

	double cumulateError(0), maxError(0);
	for (size_t w = 0; w < result.size(); ++w)
	{
		for (size_t v = 0; v < result[w].size(); ++v)
		{
			double error = utl::computeError(result[w][v], data[indices[w][v]]);
			cumulateError += error;
			maxError = (error > maxError) ? error : maxError;
		}
	}

	std::cout << "Using the mapping, the mean error is: " << utl::radians2degrees(cumulateError / nbVectors)
			  << " and the max error is: " << utl::radians2degrees(maxError) << std::endl;

	///--------------------------------------------------------------------------
	/// For comparison we try spherical Fibonacci point set quantization without 
	/// the mapping.
	///--------------------------------------------------------------------------

	cumulateError = 0;
	maxError = 0;
	for (auto & v : data)
	{
		const double error = utl::computeError(v, fib::SF(fib::inverseSF(v, nbFiboPoints), nbFiboPoints));
		cumulateError += error;
		maxError = (error > maxError) ? error : maxError;
	}
	std::cout << "Without the mapping, the mean error is: " << utl::radians2degrees(cumulateError / nbVectors) 
			  << " and the max error is: " << utl::radians2degrees(maxError) << std::endl;
	
	system("pause");
	
	
	return 0;
}