#include<algorithm>

#include "../headers/Distribution.h"
#include "../headers/MathTools.h"

namespace MathTools {

	// *********** DISTRIBUTION FUNCTIONS *********** //

	double Distribution::evaluate_CDF(const double x, const double steps) const {
		// Evaluates the CDF of the Distribution at the point x.
		//
		// This is done by numerical integration of the PDF.

		const std::vector<double> x_list = linspace(m_mean - 10 * m_stdev, x, steps);

		std::vector<double> result;
		double h;

		result.push_back(0);

		for (size_t i = 1; i < x_list.size(); i++) {
			h = x_list[i] - x_list[i - 1];
			result.push_back(result[i - 1] + 0.5 * h *
				(evaluate_PDF(x_list[i]) + evaluate_PDF(x_list[i - 1])));
		}

		return *result.end();
	}

	// *********** NORMAL DISTRIBUTION FUNCTIONS *********** //
	double NormalCDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Normal Distribution at the point x.
		//
		// See the following link for more details: https://cplusplus.com/reference/random/normal_distribution/

		return exp(-1 * pow(x - m_mean, 2) / (2 * pow(m_stdev, 2))) / (m_stdev * sqrt(2 * M_PI));
	}

	double NormalCDistribution::evaluate_inverse_CDF(const double u) const {
		// Computes the inverse CDF of the Normal distribution at the point u.
		//
		// Where the input "u" takes the values [0, 1]. Otherwise, this function returns a NaN.
		// 
		// This expression for the inverse CDF was found in Table 3.1 (Page 151) of the following book:
		// G. Fishman,Monte Carlo: Concepts, Algorithms, and Applications, Springer-Verlag, NewYork, 1996
		// Available at: https://link.springer.com/book/10.1007/978-1-4757-2553-7

		if (u < 0 || u > 1) {
			return NAN;
		}

		const double c0 = 2.515517;
		const double c1 = 0.802853;
		const double c2 = 0.010328;
		const double d1 = 1.432788;
		const double d2 = 0.189269;
		const double d3 = 0.001308;

		const double t = sqrt(-log(pow(std::min(u, 1 - u), 2)));
		
		double sign = u - 0.5;
		if (sign != 0) { sign = sign / abs(sign); }

		double result = - ( c0 + c1 * t + c2 * pow(t, 2)) / (1 + d1 * t + d2 * pow(t, 2) + d3 * pow(t,3));
		result += t;
		result *= sign * m_stdev;
		result += m_mean;

		return result;
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

	double UniformCDistribution::evaluate_inverse_CDF(const double u) const {
		// Computes the inverse CDF of the Uniform distribution at the point u.
		//
		// x = lower_lim + (upper_lim - lower_lim) * u
		//
		// Where the input "u" takes the values [0, 1]. Otherwise, this function returns a NaN.
		// 
		// This expression for the inverse CDF was found in Table 3.1 (Page 151) of the following book:
		// G. Fishman,Monte Carlo: Concepts, Algorithms, and Applications, Springer-Verlag, NewYork, 1996
		// Available at: https://link.springer.com/book/10.1007/978-1-4757-2553-7

		if (u < 0 || u > 1) {
			return NAN;
		}

		return m_lower_lim + (m_upper_lim - m_lower_lim) * u;
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

		if (x <= 0) { return 0; }

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

	// *********** GAMMA DISTRIBUTION (WIENER-ASKEY) FUNCTIONS *********** //
	void WAGammaCDistribution::set_params_stdev(const double stdev) {
		// Initialize the parameters of the Gamma Distribution using the
		// standard deviation.
		//
		// mean =  (stdev)^2
		// alpha = mean - 1
		//
		// See the following link for more details: https://doi.org/10.1137/S1064827501387826

		m_mean = stdev * stdev;
		m_stdev = stdev;

		m_alpha = m_mean - 1;
		m_alpha_shape = m_alpha + 1;
	}

	void WAGammaCDistribution::set_params_mean(const double mean) {
		// Initialize the parameters of the Gamma Distribution using the
		// standard deviation.
		//
		// stdev = (mean)^0.5
		// alpha = mean - 1
		//
		// See the following link for more details: https://doi.org/10.1137/S1064827501387826

		m_mean = mean;
		m_stdev = sqrt(mean);

		m_alpha = m_mean - 1;
		m_alpha_shape = m_alpha + 1;
	}

	double WAGammaCDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Gamma Distribution at the point x.
		//
		// See the following link for more details: https://doi.org/10.1137/S1064827501387826

		if (x <= 0) { return 0; }

		return exp(-x) * pow(x, m_alpha) / tgamma(m_alpha + 1);
	}

	double WAGammaCDistribution::sample(std::default_random_engine& generator) const {
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

		if (x < 0 || x > 1) { return 0; }

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

	// *********** BETA (WIENER-ASKEY) DISTRIBUTION FUNCTIONS *********** //
			// Constructor that uses the following equations to set the mean and
		// standard deviation
		//
		// mean = alpha / (alpha + beta)
		// stdev = sqrt(alpha * beta / (pow(alpha + beta, 2) * (alpha + beta + 1)))
		//
		// See the following link for more details: https://en.wikipedia.org/wiki/Beta_distribution
	WABetaCDistribution::WABetaCDistribution(const double alpha_new, const double beta_new) {
		m_alpha_new = alpha_new;
		m_beta_new = beta_new;
		m_alpha = m_beta_new + 1;
		m_beta = m_alpha_new + 1;

		m_mean = m_alpha / (m_alpha + m_beta);
		m_stdev = sqrt(m_alpha * m_beta / (pow(m_alpha + m_beta, 2) * (m_alpha + m_beta + 1)));

		m_mean = 2 * (m_mean - 0.5);
		m_stdev = 2 * m_stdev;
	};

	void WABetaCDistribution::set_params_mean_stdev(const double mean_new, const double stdev_new) {
		// Initialize the parameters of the Beta Distribution using the mean
		// and standard deviation. It inverts the below equations:
		//
		// mean = alpha / (alpha + beta)
		// stdev = sqrt(alpha * beta / (pow(alpha + beta, 2) * (alpha + beta + 1)))
		//
		// See the following link for more details: https://en.wikipedia.org/wiki/Beta_distribution
		//
		// This distribution has been shifted to have a mean of zero and a domain of [-1, 1].
		// In this alpha_new and beta_new are related to beta and alpha respectively.
		// See for more details: https://doi.org/10.1137/S106482750138782

		m_mean = mean_new;
		m_stdev = stdev_new;

		double mean_old = m_mean / 2. + 0.5;
		double stdev_old = m_stdev / 2.;

		m_alpha = (1 - mean_old) * pow(mean_old / stdev_old, 2) - mean_old;
		m_beta = m_alpha / mean_old * (1 - mean_old);

		if (m_alpha < 1e-3) { m_alpha = 1e-3; }
		if (m_beta < 1e-3) { m_beta = 1e-3; }

		m_alpha_new = m_beta - 1;
		m_beta_new = m_alpha - 1;
	}

	double WABetaCDistribution::evaluate_PDF(const double x) const {
		// Evaluates the PDF of the Beta Distribution at the point x.
		//
		// This distribution has been shifted to have a mean of zero and a domain of [-1, 1].
		// In this alpha_new and beta_new are related to beta and alpha respectively.
		// See for more details: https://doi.org/10.1137/S106482750138782

		if (x < -1 || x > 1) { return 0; }

		return pow(1 - x, m_alpha_new) * pow(1 + x, m_beta_new) * tgamma(m_alpha_new + m_beta_new + 2)
			/ (tgamma(m_alpha_new + 1) * tgamma(m_beta_new + 1) * pow(2, m_alpha_new + m_beta_new + 1));
	}

	double WABetaCDistribution::sample(std::default_random_engine& generator) const {
		// Takes a random sample of the Beta Distribution using a random number generator.
		//
		// See the following link for more details: https://stackoverflow.com/a/10359049
		//
		// This distribution has been shifted to have a mean of zero and a domain of [-1, 1].
		// In this alpha_new and beta_new are related to beta and alpha respectively.
		// See for more details: https://doi.org/10.1137/S106482750138782

		std::gamma_distribution<double> gammaX(m_alpha, 1);
		std::gamma_distribution<double> gammaY(m_beta, 1);

		const double X = gammaX(generator);
		const double Y = gammaY(generator);

		const double Z = X / (X + Y);

		return 2 * (Z - 0.5);
	}
}
