#include "../headers/Distribution.h"

namespace MathTools {
	double NormalDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Normal Distribution at the point x.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/

		return exp(-1 * pow(x - m_mean, 2) / (2 * pow(m_stdev, 2))) / (m_stdev * pow(2 * M_PI, 0.5));
	}

	double NormalDistribution::sample(std::default_random_engine& generator) const {
		// Takes a random sample of the normal distribution using a random number generator
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/

		std::normal_distribution<double> norm(m_mean, m_stdev);

		return norm(generator);
	}
}
