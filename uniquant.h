// This file is part of the reference implementation for the paper
//    Fast Lossy Compression of 3D Unit Vector Sets
//    Sylvain Rousseau and Tamy Boubekeur
//    SIGGRAPH Asia 2017 Technical Briefs (SA '17)
//    DOI: https://doi.org/10.1145/3145749.3149436
//
// All rights reserved. Use of this source code is governed by a 
// MIT license that can be found in the LICENSE file.

#pragma once

#include "externals/glm/glm.hpp"

//#define __BMI2__
#include "externals/libmorton/include/morton.h"

//linux compatibility
#undef M_PI
#undef M_TWO_PI

/// 
///	The namespace uvc contain the code necessary to the compression algorithm
///
namespace uvc
{
	const double M_PI	  = 3.1415926535897932384626433832795028841971693993751058209749445923;
	const double M_TWO_PI = 6.2831853071795864769252867665590057683943387987502116419498891846;
	
	/// convert cartesian coordinates of an unit vector to spherical coordinates
	inline glm::fvec2 cartesianToSpherical(const glm::fvec3 & vector)
	{
		return glm::fvec2(std::acos(vector[2]), glm::atan(vector[1], vector[0]));
	}

	/// convert spherical coordinates to an cartesian unit vector
	inline glm::fvec3 sphericalToCartesian(const glm::fvec2 & unitVector)
	{
		double sinx = -sin(unitVector[0]);
		return glm::fvec3(sinx * cos(unitVector[1]), sinx * sin(unitVector[1]), cos(unitVector[0]));
	}

	/// calculate the morton code of the discretized spherical coordinates of the unit vector $vector
	inline uint32_t mortonCodeUnitVector(const glm::fvec3 & vector)
	{
		glm::fvec2 v = cartesianToSpherical(vector);
		uint16_t x = uint16_t(0xffff * (v[0] / M_PI));
		uint16_t y = uint16_t(0xffff * ((v[1] + M_PI) / M_TWO_PI));
		return morton2D_32_encode(x, y);
	}

	/// group unit vectors contained in $data according to the $level most significant bits of the morton code of their discretized spherical coordinates. 
	/// As we are working on 32 bits morton code, $level needs to be in ]1, 32[. The output is $groupedData containing for each group the vectors and $indices
	/// An array making possible to retreive the position of the vector in $data.	 
	inline void group(const std::vector<glm::fvec3> & data, 
		std::vector<std::vector<size_t>> & indices, 
		std::vector<std::vector<glm::fvec3>> & groupedData, 
		const unsigned level)
	{
		/// we are working on 32 bits
		assert(level < 32); 
		assert(level > 1);

		indices.clear();
		groupedData.clear();
		
		/// allocate memory
		uint32_t nbGroups = uint32_t(1 << level);
		size_t meanSize = data.size() / nbGroups;

		groupedData.resize(nbGroups);
		indices.resize(nbGroups);
		
		for (unsigned i = 0; i < indices.size(); ++i)
		{
			indices[i].reserve(meanSize);
		}

		/// add vectors to the correct window, register to the indice array
		for (unsigned i = 0; i < data.size(); ++i)
		{
			unsigned id = (mortonCodeUnitVector(data[i]) >> (32 - level));
			groupedData[id].push_back(data[i]);
			indices[id].push_back(i);
		}
	}

	/// Mapping from a point on the surface of the unit sphere to a point on the spherical cap centered in $axis 
	/// whose projection on $axis is $ratio. The mapping is such that a point set with an uniform distribution 
	/// on the unit sphere will be mapped to a point set with a uniform distribution on the surface of the spherical cap (cf 2.3).
	inline glm::fvec3 uniformMapping(const glm::fvec3 & v, const glm::fvec3 & axis, const float ratio)
	{
		glm::dvec3 v2(v);
		glm::dvec3 axis2(axis);
		v2 = glm::normalize(v2);
		axis2 = glm::normalize(axis2);

		double K = 1.0 / ratio;
		double c = glm::dot(v2, axis2);
		
		/// numerical instabilities
		if (c > 1.0)
			c = 1.0;
		if (c < -1.0)
			c = -1.0;
		glm::dvec3 p1 = glm::normalize(v2 - c * axis2);
		double delta = (1.0 - ((1.0 - c) / K));
		return  p1 * std::sqrt(1.0 - delta * delta) + delta * axis2;
	}

	/// Mapping from a point on the spherical cap centered in $axis whose projection on $axis is $ratio.
	/// to a point on the surface of the unit sphere. The mapping is such that a point set with an uniform distribution 
	/// on the spherical cap will be mapped to a point set with a uniform distribution on the surface of the unit sphere (cf 2.3).
	///
	inline glm::fvec3 inverseUniformMapping(const glm::fvec3 & v, const glm::fvec3 & axis, const float ratio)
	{
		glm::dvec3 v2(v);
		glm::dvec3 axis2(axis);
		v2 = glm::normalize(v2);
		axis2 = glm::normalize(axis2);

		double K = ratio;
		double c = glm::dot(v2, axis2);
		///numerical instabilities, we need to check $c
		if (c > 1.0)
			c = 1.0;
		if (c < -1.0)
			c = -1.0;
		glm::dvec3 p1 = glm::normalize(v2 - c * axis2);
		double delta = (1.0 - ((1.0 - c) / K));
		return  p1 * std::sqrt(1.0 - delta * delta) + delta * axis2;
	}

	/// create the mask with $nbBits 1 and right zero padding
	constexpr uint32_t NMostSignificantBitsMask(const unsigned nbBits)
	{
		return ~(0xffffffff >> nbBits);
	}

	/// Compute a unit vector whose morton code is $mortonCode
	inline glm::fvec3 unitVectorFromMortonCode(const uint32_t mortonCode)
	{
		uint_fast16_t x, y;
		morton2D_32_decode(mortonCode, x, y);
		glm::fvec2 sphericalCoordinates((double(x) * M_PI) / double(0xffff), ((double(y) - M_PI) * (2.f * M_PI)) / double(0xffff));
		return glm::normalize(sphericalToCartesian(sphericalCoordinates));
	}

	/// compute the average vector of the group containing v with the mask of the group being $mask
	inline glm::fvec3 computeAverage(const glm::fvec3 & v, const uint32_t mask)
	{
		uint32_t mc = (mortonCodeUnitVector(v) & mask) + ((~mask) / 4);
		return unitVectorFromMortonCode(mc);
	}

	/// Compute the average vector from the group id and the mask
	inline glm::fvec3 computeAverage(const uint32_t groupID, const uint32_t mask, const uint32_t level)
	{
		uint32_t mc = (groupID << (32 - level)) + ((~mask) / 4);
		return unitVectorFromMortonCode(mc);
	}

	/// compute the ratio (K in Eq. 19) for the mapping of a group containing vector 
	inline float computeRatio(const glm::fvec3 & vector, const uint32_t mask)
	{
		uint32_t key = mortonCodeUnitVector(vector) & mask;
		float mindot = glm::dot(vector, unitVectorFromMortonCode(key));
		float dot = glm::dot(vector, unitVectorFromMortonCode(key | (~mask)));
		mindot = (dot < mindot) ? dot : mindot;
		return ((1.f - mindot) / 2.f) + 0.005f; /// 0.005f is the epsilone (cf 2.4)
	}
}
