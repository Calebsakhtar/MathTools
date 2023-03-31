#ifndef CDISTRIBUTION_H
#define CDISTRIBUTION_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <string>
#include <random>

namespace MathTools {

	// *********** GENERIC DISTRIBUTION (PARENT CLASS) *********** //
	class CDistribution {
	protected:
		std::string m_name;
		double m_mean;
		double m_stdev;

	public:
		CDistribution() {};

		CDistribution(const double mean, const double stdev) {
			m_mean = mean;
			m_stdev = stdev;
		}

		double evaluate_PDF(const double x) const;

		double sample(std::default_random_engine& generator) const;
	};

	// *********** NORMAL DISTRIBUTION *********** //
	class NormalCDistribution: public CDistribution {
		// Use the same constructors as CDistribution
		using CDistribution::CDistribution;

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
	class UniformCDistribution : public CDistribution {
		double m_lower_lim; // Lower limit for the distribution
		double m_upper_lim; // Upper limit for the distribution

	public:
		UniformCDistribution() {};

		// Constructor that uses the following equations to set the mean and
		// standard deviation
		//
		// mean = (lower_lim + upper_lim)/2
		// stdev = (upper_lim - lower_lim)/sqrt(12)
		//
		// See the following link for more details: https://www.statology.org/uniform-distribution-r/#:~:text=The%20uniform%20distribution%20has%20the,distribution%20is%20%CF%83%20%3D%20%E2%88%9A%CF%83
		UniformCDistribution(const double lower_lim, const double upper_lim) {
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

	// *********** GAMMA DISTRIBUTION *********** //
	class GammaCDistribution : public CDistribution {
		double m_alpha_shape; // Shape parameter
		double m_beta_scale; // Scale parameter

	public:
		GammaCDistribution() {};

		// Constructor that uses the following equations to set the mean and
		// standard deviation
		//
		// mean = alpha/beta
		// stdev = sqrt(alpha)/beta
		//
		// See the following link for more details: https://math.stackexchange.com/questions/1810257/gamma-functions-mean-and-standard-deviation-through-shape-and-rate#:~:text=A%20gamma%20distribution%20has%20a,%5D%3D%E2%88%9Aa%2Fb.
		GammaCDistribution(const double alpha_shape, const double beta_scale) {
			m_alpha_shape = alpha_shape;
			m_beta_scale = beta_scale;

			m_mean = alpha_shape / beta_scale;
			m_stdev = sqrt(alpha_shape) / beta_scale;
		};

		// Initialize the parameters of the Gamma Distribution using the mean
		// and standard deviation.
		//
		// alpha = (mean/stdev)^2
		// beta = mean/stdev^2
		//
		// See the following link for more details: https://math.stackexchange.com/questions/1810257/gamma-functions-mean-and-standard-deviation-through-shape-and-rate#:~:text=A%20gamma%20distribution%20has%20a,%5D%3D%E2%88%9Aa%2Fb.
		void set_params_mean_stdev(const double mean, const double stdev);

		// Evaluates the PDF of the Gamma Distribution at the point x.
		//
		// See the following link for more details: https://en.cppreference.com/w/cpp/numeric/random/gamma_distribution
		virtual double evaluate_PDF(const double x) const; // Virtual overrides parent class method

		// Takes a random sample of the Gamma Distribution using a random number generator.
		//
		// See the following link for more details: https://en.cppreference.com/w/cpp/numeric/random/gamma_distribution
		virtual double sample(std::default_random_engine& generator) const;
	};

	// *********** BETA DISTRIBUTION *********** //
	class BetaCDistribution : public CDistribution {
		double m_alpha;
		double m_beta;

	public:
		BetaCDistribution() {};

		// Constructor that uses the following equations to set the mean and
		// standard deviation
		//
		// mean = alpha / (alpha + beta)
		// stdev = sqrt(alpha * beta / (pow(alpha + beta, 2) * (alpha + beta + 1)))
		//
		// See the following link for more details: https://en.wikipedia.org/wiki/Beta_distribution
		BetaCDistribution(const double alpha, const double beta) {
			m_alpha = alpha;
			m_beta = beta;

			m_mean = alpha / (alpha + beta);
			m_stdev = sqrt(alpha * beta /( pow(alpha + beta,2) * (alpha + beta + 1)));
		};

		// Initialize the parameters of the Beta Distribution using the mean
		// and standard deviation. It inverts the below equations:
		//
		// mean = alpha / (alpha + beta)
		// stdev = sqrt(alpha * beta / (pow(alpha + beta, 2) * (alpha + beta + 1)))
		//
		// See the following link for more details: https://en.wikipedia.org/wiki/Beta_distribution
		void set_params_mean_stdev(const double mean, const double stdev);

		// Evaluates the PDF of the Beta Distribution at the point x.
		//
		// See the following link for more details: 
		virtual double evaluate_PDF(const double x) const; // Virtual overrides parent class method

		// Takes a random sample of the Beta Distribution using a random number generator.
		//
		// See the following link for more details: https://stackoverflow.com/a/10359049
		virtual double sample(std::default_random_engine& generator) const;
	};
}

#endif
