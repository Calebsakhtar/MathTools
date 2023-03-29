#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <string>
#include <random>

namespace MathTools {

	// *********** GENERIC DISTRIBUTION (PARENT CLASS) *********** //
	class Distribution {
	protected:
		std::string m_name;
		double m_mean;
		double m_stdev;

	public:
		Distribution() {};

		Distribution(const double mean, const double stdev) {
			m_mean = mean;
			m_stdev = stdev;
		}

		double evaluate_PDF(const double x) const;

		double sample(std::default_random_engine& generator) const;
	};

	// *********** NORMAL DISTRIBUTION *********** //
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

	// *********** UNIFORM DISTRIBUTION *********** //
	class UniformDistribution : public Distribution {
		double m_lower_lim; // Lower limit for the distribution
		double m_upper_lim; // Upper limit for the distribution

	public:
		UniformDistribution() {};

		// Constructor that uses the following equations to set the mean and
		// standard deviation
		//
		// mean = (lower_lim + upper_lim)/2
		// stdev = (upper_lim - lower_lim)/sqrt(12)
		//
		// See the following link for more details: https://www.statology.org/uniform-distribution-r/#:~:text=The%20uniform%20distribution%20has%20the,distribution%20is%20%CF%83%20%3D%20%E2%88%9A%CF%83
		UniformDistribution(const double lower_lim, const double upper_lim) {
			m_lower_lim = lower_lim;
			m_upper_lim = upper_lim;

			m_mean = (lower_lim + upper_lim) * 0.5;
			m_stdev = (upper_lim - lower_lim) / sqrt(12.0);
		};

		// Initialize the parameters of the uniform distribution using the mean
		// and standard deviation.
		//
		// mean = (lower_lim + upper_lim)/2
		// stdev = (upper_lim - lower_lim)/sqrt(12)
		//
		// See the following link for more details: https://www.statology.org/uniform-distribution-r/#:~:text=The%20uniform%20distribution%20has%20the,distribution%20is%20%CF%83%20%3D%20%E2%88%9A%CF%83
		void set_params_mean_stdev(const double mean, const double stdev);

		// Evaluates the PDF of the Uniform Distribution at the point x.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/uniform_real_distribution/
		virtual double evaluate_PDF(const double x) const; // Virtual overrides parent class method

		// Takes a random sample of the Uniform Distribution using a random number generator.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/uniform_real_distribution/
		virtual double sample(std::default_random_engine& generator) const;
	};

}

#endif
