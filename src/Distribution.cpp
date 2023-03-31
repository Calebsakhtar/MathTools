#include "../headers/Distribution.h"

namespace MathTools {

	// *********** NORMAL DISTRIBUTION FUNCTIONS *********** //
	double NormalCDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Normal Distribution at the point x.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/

		return exp(-1 * pow(x - m_mean, 2) / (2 * pow(m_stdev, 2))) / (m_stdev * sqrt(2 * M_PI));
	}

	double NormalCDistribution::sample(std::default_random_engine& generator) const {
		// Takes a random sample of the normal distribution using a random number generator
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/

		std::normal_distribution<double> norm(m_mean, m_stdev);

		return norm(generator);
	}

	// *********** UNIFORM DISTRIBUTION FUNCTIONS *********** //
	void UniformCDistribution::set_params_mean_stdev(const double mean, const double stdev) {
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

	double UniformCDistribution::evaluate_PDF(const double x) const {
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

	double UniformCDistribution::sample(std::default_random_engine& generator) const {
		// Takes a random sample of the Uniform Distribution using a random number generator.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/uniform_real_distribution/

		std::uniform_real_distribution<double> unif(m_lower_lim, m_upper_lim);

		return unif(generator);
	}

	// *********** GAMMA DISTRIBUTION FUNCTIONS *********** //
	void GammaCDistribution::set_params_mean_stdev(const double mean, const double stdev) {
		// Initialize the parameters of the Gamma Distribution using the mean
		// and standard deviation.
		//
		// alpha = (mean/stdev)^2
		// beta = mean/stdev^2
		//
		// See the following link for more details: https://math.stackexchange.com/questions/1810257/gamma-functions-mean-and-standard-deviation-through-shape-and-rate#:~:text=A%20gamma%20distribution%20has%20a,%5D%3D%E2%88%9Aa%2Fb.

		m_mean = mean;
		m_stdev = stdev;

		m_alpha_shape = pow(mean / stdev, 2);
		m_beta_scale = mean / pow(stdev, 2);
	}

	double GammaCDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Gamma Distribution at the point x.
		//
		// See the following link for more details: https://en.cppreference.com/w/cpp/numeric/random/gamma_distribution

		return exp(-x / m_beta_scale) * pow(x, m_alpha_shape - 1) / 
			(pow(m_beta_scale, m_alpha_shape) * tgamma(m_alpha_shape));
	}

	double GammaCDistribution::sample(std::default_random_engine& generator) const {
		// Takes a random sample of the Gamma Distribution using a random number generator.
		//
		// See the following link for more details: https://en.cppreference.com/w/cpp/numeric/random/gamma_distribution

		std::gamma_distribution<double> gamma(m_alpha_shape, m_beta_scale);

		return gamma(generator);
	}

	// *********** BETA DISTRIBUTION FUNCTIONS *********** //
	void BetaCDistribution::set_params_mean_stdev(const double mean, const double stdev){
		// Initialize the parameters of the Beta Distribution using the mean
		// and standard deviation. It inverts the below equations:
		//
		// mean = alpha / (alpha + beta)
		// stdev = sqrt(alpha * beta / (pow(alpha + beta, 2) * (alpha + beta + 1)))
		//
		// See the following link for more details: https://en.wikipedia.org/wiki/Beta_distribution

		m_alpha = (1 - mean) * pow(mean / stdev, 2) - mean;
		m_beta = m_alpha / mean * (1 - mean);
	}

	double BetaCDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Beta Distribution at the point x.
		//
		// See the following link for more details: https://en.wikipedia.org/wiki/Beta_distribution

		return pow(x, m_alpha - 1) * pow(1 - x, m_beta - 1)	* tgamma(m_alpha + m_beta)
			/ (tgamma(m_alpha) * tgamma(m_beta));
	}

	double BetaCDistribution::sample(std::default_random_engine& generator) const {
		// Takes a random sample of the Beta Distribution using a random number generator.
		//
		// See the following link for more details: https://stackoverflow.com/a/10359049

		std::gamma_distribution<double> gammaX(m_alpha, 1);
		std::gamma_distribution<double> gammaY(m_beta, 1);

		const double X = gammaX(generator);
		const double Y = gammaY(generator);

		return X / (X + Y);
	}
}
