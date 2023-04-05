#ifndef ORTH_POLY_H
#define ORTH_POLY_H

#include <cmath>
#include <vector>

namespace MathTools {

	class OrthPoly {
	protected:
		unsigned int m_n; // Polynomial order

	public:
		OrthPoly() {};

		OrthPoly(const unsigned int n) {
			m_n = n;
		};

		double evaluate(const double x) const;

		unsigned int get_order() const{
			return m_n;
		}

		void set_order(const unsigned int n){
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
}

#endif // !ORTH_POLY_H
