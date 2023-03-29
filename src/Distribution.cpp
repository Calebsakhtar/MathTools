#include "../headers/Distribution.h"

namespace MathTools {

	// *********** NORMAL DISTRIBUTION FUNCTIONS *********** //
	double NormalDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Normal Distribution at the point x.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/

		return exp(-1 * pow(x - m_mean, 2) / (2 * pow(m_stdev, 2))) / (m_stdev * sqrt(2 * M_PI));
	}

	double NormalDistribution::sample(std::default_random_engine& generator) const {
		// Takes a random sample of the normal distribution using a random number generator
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/

		std::normal_distribution<double> norm(m_mean, m_stdev);

		return norm(generator);
	}

	// *********** UNIFORM DISTRIBUTION FUNCTIONS *********** //
	void UniformDistribution::set_params_mean_stdev(const double mean, const double stdev) {
		// Initialize the parameters of the uniform distribution using the mean
		// and standard deviation.
		//
		// mean = (lower_lim + upper_lim)/2
		// stdev = (upper_lim - lower_lim)/sqrt(12)
		//
		// See the following link for more details: https://www.statology.org/uniform-distribution-r/#:~:text=The%20uniform%20distribution%20has%20the,distribution%20is%20%CF%83%20%3D%20%E2%88%9A%CF%83
		
		m_mean = mean;
		m_stdev = stdev;

		const double delta = sqrt(3) * stdev;
		m_lower_lim = mean - delta;
		m_upper_lim = mean + delta;
	}

	double UniformDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Uniform Distribution at the point x.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/uniform_real_distribution/

		if ((x < m_lower_lim) || (x >= m_upper_lim)) {
			return 0.0;
		}
		else {
			return 1.0 / (m_upper_lim - m_lower_lim);
		}
	}

	double UniformDistribution::sample(std::default_random_engine& generator) const {
		// Takes a random sample of the Uniform Distribution using a random number generator.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/uniform_real_distribution/

		std::uniform_real_distribution<double> unif(m_lower_lim, m_upper_lim);

		return unif(generator);
	}
}
