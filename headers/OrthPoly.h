#ifndef ORTH_POLY_H
#define ORTH_POLY_H

#include <cmath>
#include <vector>

namespace MathTools {

	class OrthPoly {
	protected:
		size_t m_n; // Polynomial order
		std::vector<double> m_coeffs = {};

	public:
		OrthPoly() {};

		OrthPoly(const size_t n) {
			m_n = n;
		};

		virtual double evaluate(const double x) const = 0;

		size_t get_order() const{
			return m_n;
		}

		void set_order(const size_t n){
			m_n = n;
		};
		
	};

	// *********** LEGENDRE POLYNOMIALS *********** //
	class LegendrePoly : public OrthPoly {
	public:
		LegendrePoly() {};

		// Use the summation expression for the Legendre polynomials
		// to store the relevant coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas
		LegendrePoly(const size_t n);

		// Evaluate the summation form of the Legendre polynomial 
		// using the stored coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas
		double evaluate(const double x) const;
	};

	// *********** HERMITE POLYNOMIALS *********** //
	class HermitePoly : public OrthPoly {
	public:
		HermitePoly() {};

		// Use the summation expression for the Hermite polynomials
		// to store the relevant coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Hermite_polynomials#Explicit_expression
		HermitePoly(const size_t n);

		double evaluate(const double x) const;
	};

	// *********** LAGUERRE POLYNOMIALS *********** //
	class LaguerrePoly : public OrthPoly {
	public:
		LaguerrePoly() {};

		// Use the summation expression for the Laguerre polynomials
		// to store the relevant coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Laguerre_polynomials#Recursive_definition,_closed_form,_and_generating_function
		LaguerrePoly(const size_t n);

		// Evaluate the summation form of the Laguerre polynomial 
		// using the stored coefficients.
		//
		// See: https://en.wikipedia.org/wiki/Laguerre_polynomials#Recursive_definition,_closed_form,_and_generating_function
		double evaluate(const double x) const;
	};

	// *********** JACOBI POLYNOMIALS *********** //
	class JacobiPoly : public OrthPoly {
		double m_alpha = 0; // parameter alpha
		double m_beta = 0; // parameter beta

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
		double evaluate(const double x) const;

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
