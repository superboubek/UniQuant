// This file is part of the reference implementation for the paper
//    Fast Lossy Compression of 3D Unit Vector Sets
//    Sylvain Rousseau and Tamy Boubekeur
//    SIGGRAPH Asia 2017 Technical Briefs (SA '17)
//    DOI: https://doi.org/10.1145/3145749.3149436
//
// All rights reserved. Use of this source code is governed by a 
// MIT license that can be found in the LICENSE file.


/// implementation of spherical fibonacci and inverse spherical fibonacci inspired by 
/// https://www.shadertoy.com/view/lllXz4 witch using the method described in the paper: 
/// Spherical Fibonacci Mapping
/// Benjamin Keinert, Matthias Innmann, Michael Sänger and Marc Stamminger
/// Proc. of SIGGRAPH Asia 2015
/// The original source code was realized by Inigo Quilez and was released under the MIT license

#pragma once


#include "externals/glm/glm.hpp"

namespace fib
{
	const double M_PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
	const double M_TWO_PI = 6.2831853071795864769252867665590057683943387987502116419498891846;
	const double PHI = 1.6180339887498948482045868343656381177203091798057628621354486227;

	/// get the coordinates of point $id of the spherical fibonacci point set containing $nbPoints points
	glm::dvec3 SF(const uint32_t id, const uint32_t nbPoints)
	{
		double m = 1.0 - 1.0 / nbPoints;
		double phi = M_TWO_PI * glm::fract(id * PHI);
		double cosTheta = m - 2.0 * double(id) / nbPoints;
		double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
		return glm::dvec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
	}

	/// get the id of the closest point to $p$ in the spherical fibonacci point set containing $n$ points
	/// due to numerical instabilities, keep n under 2^24
	uint32_t inverseSF(const glm::dvec3 & p, const unsigned n)
	{
		double m = 1.0 - 1.0 / n;
		double phi = glm::min(glm::atan(p.y, p.x), M_PI), cosTheta = p.z;

		double k = glm::max(2.0, floor(log(n * M_PI * sqrt(5.0) * (1.0 - cosTheta*cosTheta)) / log(PHI + 1.0)));
		double Fk = pow(PHI, k) / sqrt(5.0);
		glm::dvec2  F = glm::vec2(round(Fk), round(Fk * PHI));

		glm::vec2 ka = 2.0 * F / double(n);
		glm::vec2 kb = M_TWO_PI *(glm::fract((F + 1.0) * PHI) - (PHI - 1.0));

		glm::mat2 iB = glm::mat2(ka.y, -ka.x,
			kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);

		glm::dvec2 c = floor(iB * glm::vec2(phi, cosTheta - m));
		double d = 8.0;
		uint32_t j = 0;
		for (int s = 0; s < 4; s++)
		{
			glm::dvec2 uv = glm::dvec2(double(s - 2 * (s / 2)), double(s / 2));
			double i = glm::dot(F, uv + c);

			double phi = M_TWO_PI * glm::fract(i * PHI);
			double cosTheta = m - 2.0 * i / n;
			double sinTheta = sqrt(1.0 - cosTheta*cosTheta);

			glm::dvec3 q = glm::vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
			double squaredDistance = glm::dot(q - p, q - p);
			if (squaredDistance < d)
			{
				d = squaredDistance;
				j = uint32_t(i);
			}
		}
		return j;
	}

}