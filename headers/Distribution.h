#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <string>
#include <random>

namespace MathTools {

	// GENERIC DISTRIBUTION (PARENT CLASS)
	class Distribution {
	protected:
		std::string m_name;
		double m_mean;
		double m_stdev;

	public:

		Distribution(const double mean, const double stdev) {
			m_mean = mean;
			m_stdev = stdev;
		}

		double evaluate_PDF(const double x) const;

		double sample(std::default_random_engine& generator) const;
	};

	// NORMAL DISTRIBUTION
	class NormalDistribution: public Distribution {
		// Use the same constructors as Distribution
		using Distribution::Distribution;

	public:
		// Evaluates the PDF of the Normal Distribution at the point x.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/
		virtual double evaluate_PDF(const double x) const; // Virtual overrides parent class method

		// Takes a random sample of the normal distribution using a random number generator
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/
		virtual double sample(std::default_random_engine& generator) const;
	};


}

#endif
