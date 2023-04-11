#include "../headers/OrthPoly.h"
#include "../headers/MathTools.h"

namespace MathTools {

	JacobiPoly::JacobiPoly(const size_t n, const double alpha,
		// Use the summation expression for the Jacobi polynomials
		// to store the relevant coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Jacobi_polynomials#Alternate_expression_for_real_argument

		const double beta) {
		m_n = n;
		m_alpha = alpha;
		m_beta = beta;
 
		double coeff = 0;
		double s_double = 0;
		const double n_double = static_cast<double>(n);

		// Loop from 0 to n, both inclusive
		for (size_t s = 0; s < n + 1; s++){
			s_double = static_cast<double>(s);
			coeff = nChoosek_gamma(n_double + alpha, n_double - s_double);
			coeff /= pow(2, n_double);
			coeff *= nChoosek_gamma(n_double + beta, s_double);
			m_coeffs.push_back(coeff);
		}
	}

	double JacobiPoly::evaluate(const double x) const {
		// Evaluate the summation form of the Jacobi polynomial 
		// using the stored coefficient.
		//
		// See: https://en.wikipedia.org/wiki/Jacobi_polynomials#Alternate_expression_for_real_argument

		double result = 0;
		double s_double = 0;
		const double n_double = static_cast<double>(m_n);

		for (size_t s = 0; s < m_coeffs.size(); s++) {
			s_double = static_cast<double>(s);
			result += m_coeffs[s] * pow(x - 1, s_double)
				* pow(x + 1, n_double - s_double);
		}

		return result;
	}
}
