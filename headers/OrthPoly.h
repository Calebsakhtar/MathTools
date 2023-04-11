#ifndef ORTH_POLY_H
#define ORTH_POLY_H

#include <cmath>
#include <vector>

namespace MathTools {

	class OrthPoly {
	protected:
		size_t m_n; // Polynomial order

	public:
		OrthPoly() {};

		OrthPoly(const size_t n) {
			m_n = n;
		};

		double evaluate(const double x) const;

		size_t get_order() const{
			return m_n;
		}

		void set_order(const size_t n){
			m_n = n;
		};
		
	};

	// *********** LEGENDRE POLYNOMIALS *********** //
	class LegendrePoly: public OrthPoly {
		// Use the same constructors as OrthPoly
		using OrthPoly::OrthPoly;
	public:
		virtual double evaluate(const double x) const{
			// See documentation: https://en.cppreference.com/w/cpp/numeric/special_functions/legendre			
			return std::legendre(m_n, x);
		}
	};

	// *********** HERMITE POLYNOMIALS *********** //
	class HermitePoly : public OrthPoly {
		// Use the same constructors as OrthPoly
		using OrthPoly::OrthPoly;
	public:
		virtual double evaluate(const double x) const{
			// See documentation: https://en.cppreference.com/w/cpp/numeric/special_functions/hermite
			return std::hermite(m_n, x);
		}
	};

	// *********** LAGUERRE POLYNOMIALS *********** //
	class LaguerrePoly : public OrthPoly {
		// Use the same constructors as OrthPoly
		using OrthPoly::OrthPoly;
	public:
		virtual double evaluate(const double x) const{
			// See documentation: https://en.cppreference.com/w/cpp/numeric/special_functions/laguerre
			return std::laguerre(m_n, x);
		}
	};

	// *********** JACOBI POLYNOMIALS *********** //
	class JacobiPoly : public OrthPoly {
		double m_alpha = 0; // parameter alpha
		double m_beta = 0; // parameter beta
		std::vector<double> m_coeffs = {};

	public:
		JacobiPoly() {};

		// Use the summation expression for the Jacobi polynomials
		// to store the relevant coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Jacobi_polynomials#Alternate_expression_for_real_argument
		JacobiPoly(const size_t n, const double alpha, const double beta);

		// Evaluate the summation form of the Jacobi polynomial 
		// using the stored coefficient.
		//
		// See: https://en.wikipedia.org/wiki/Jacobi_polynomials#Alternate_expression_for_real_argument
		virtual double evaluate(const double x) const;

		double get_alpha() const {
			return m_alpha;
		}

		void set_alpha(const double alpha) {
			m_alpha = alpha;
		};

		double get_beta() const {
			return m_beta;
		}

		void set_beta(const double beta) {
			m_alpha = beta;
		};

	};
}

#endif // !ORTH_POLY_H
