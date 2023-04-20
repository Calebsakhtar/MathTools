#include "../headers/OrthPoly.h"
#include "../headers/MathTools.h"

namespace MathTools {

	// *********** LEGENDRE POLYNOMIALS *********** //
	LegendrePoly::LegendrePoly(const size_t n) {
		// Use the summation expression for the Legendre polynomials
		// to store the relevant coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas

		m_n = n;
		const size_t k_max = n / 2;
		const double n_double = static_cast<double>(n);
		
		double power_2 = pow(2, -1 * n_double);
		double k_double = 0;
		double coeff = 0;

		// Loop from 0 to k_max, both inclusive
		for (size_t k = 0; k < k_max + 1; k++) {
			k_double = static_cast<double>(k);
			coeff = nChoosek_gamma(n_double, k_double);
			coeff *= power_2;
			coeff *= nChoosek_gamma(2 * n_double - 2 * k_double, n_double);

			// If k is odd, make the coefficient negative
			if (k % 2 == 1) { coeff *= -1; }

			m_coeffs.push_back(coeff);
		}
	}

	double LegendrePoly::evaluate(const double x) const{
		// Evaluate the summation form of the Legendre polynomial 
		// using the stored coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas

		double result = 0;
		double k_double = 0;
		const double n_double = static_cast<double>(m_n);

		for (size_t k = 0; k < m_coeffs.size(); k++) {
			k_double = static_cast<double>(k);
			result += m_coeffs[k] * pow(x, n_double - 2 * k_double);
		}

		return result;
	}

	// *********** HERMITE POLYNOMIALS *********** //
	HermitePoly::HermitePoly(const size_t n) {
		// Use the summation expression for the Hermite polynomials
		// to store the relevant coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Hermite_polynomials#Explicit_expression

		m_n = n;
		const size_t k_max = n / 2;
		const double n_double = static_cast<double>(n);

		double k_double = 0;
		double coeff = 0;

		bool finished = false;
		bool p_finished = false;
		bool q_finished = false;
		bool r_finished = false;

		// Used to compute the triple factorial product
		size_t p = 0;
		size_t q = 0;
		size_t r = 0;

		// Loop from 0 to k_max, both inclusive
		for (size_t k = 0; k < k_max + 1; k++) {
			k_double = static_cast<double>(k);

			coeff = pow(2, -1 * k_double);
			
			// If k is odd, make the coefficient negative
			if (k % 2 == 1) { coeff *= -1; }

			p = n;
			q = k;
			r = n - 2 * k;

			// Reset all the booleans
			finished = false;
			p_finished = false;
			q_finished = false;
			r_finished = false;

			// Compute the triple factorial product
			while (!finished) {
				if (p > 1) {
					coeff *= static_cast<double>(p);
					p--;
				}
				else {
					p_finished = true;
				}

				if (q > 1) {
					coeff /= static_cast<double>(q);
					q--;
				}
				else {
					q_finished = true;
				}

				if (r > 1) {
					coeff /= static_cast<double>(r);
					r--;
				}
				else {
					r_finished = true;
				}

				finished = p_finished && q_finished && r_finished;
			}

			m_coeffs.push_back(coeff);
		}
	}
	
	double HermitePoly::evaluate(const double x) const {
		// Evaluate the summation form of the Hermite polynomial 
		// using the stored coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Hermite_polynomials#Explicit_expression

		double result = 0;
		double k_double = 0;
		const double n_double = static_cast<double>(m_n);

		for (size_t k = 0; k < m_coeffs.size(); k++) {
			k_double = static_cast<double>(k);
			result += m_coeffs[k] * pow(x, n_double - 2 * k_double);
		}

		return result;
	}

	// *********** JACOBI POLYNOMIALS *********** //
	JacobiPoly::JacobiPoly(const size_t n, const double alpha,
		const double beta) {
		// Use the summation expression for the Jacobi polynomials
		// to store the relevant coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Jacobi_polynomials#Alternate_expression_for_real_argument
		
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
		// using the stored coefficients.
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
