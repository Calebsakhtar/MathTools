#include <math.h>
#include <fstream>
#include "../headers/MathTools.h"

namespace MathTools {

	std::vector<double> linspace(const double start_val, const double end_val,
		const size_t num_pts) {
		// Returns a std::vector containing num_points evenly spaced between start_val and end_val.
		//
		// num_points is inclusive of start_val and end_val
		std::vector<double> results = { start_val };
		const double increment = (end_val - start_val) / static_cast<double>(num_pts - 1);

		for (size_t i = 1; i < num_pts; i++) {
			results.push_back(results[i - 1] + increment);
		}

		return results;
	}

	bool is_same_double(const double A, const double B, const double tol) {
		// Given two doubles, this function returns whether they are within tol of each other.
		//
		// The input tol is optional, and has a default value of 1e-6.

		return abs(A - B) < tol;
	}

	void print_scalar_list(const std::vector<double>& scalar_list,
		// Prints a list of scalars into a ".csv" file

		const std::string& filename) {
		std::string fileloc = "../Outputs/";
		fileloc += filename + ".csv";

		// Create and open a text file to write into
		std::ofstream Opfile(fileloc);

		// Loop through the times vector and print its elements
		for (size_t i = 0; i < scalar_list.size(); i++) {
			Opfile << std::to_string(i) << "," << std::to_string(scalar_list[i]) << "\n";
		}

		// Close the output file
		Opfile.close();
	}

	void print_vector_list(const std::vector<Eigen::VectorXd>& vect_list,
		const std::string& filename) {
		// Prints a list of vectors into a ".csv" file

		// Create the File Location String
		std::string fileloc = "../Outputs/";
		fileloc += filename + ".csv";

		// Create and open a text file to write into
		std::ofstream Opfile(fileloc);

		// Loop through each solution vector and print its elements
		for (size_t i = 0; i < vect_list.size(); i++) {
			Opfile << std::to_string(i);
			const Eigen::VectorXd current_vect = vect_list[i];

			for (size_t j = 0; j < vect_list[0].size(); j++) {
				Opfile << "," << std::to_string(current_vect[j]);
			}

			Opfile << "\n";
		}

		// Close the output file
		Opfile.close();
	}

	template <typename int_type>
	int_type factorial(const int_type n) {
		// Returns the factorial of an integer-like input.
		// The type of the output matches that of the input.

		int_type soln = n;

		for (int_type i = n - 1; i > 1; i--){
			soln *= i;
		}

		return soln;
	}

	double factorial_gamma(const double n) {
		// Returns the value of the gamma function at the input plus one.
		// For positive integers, this is equivalent to computing the factorial.
		//
		// See for more details: https://www.cantorsparadise.com/the-beautiful-gamma-function-and-the-genius-who-discovered-it-8778437565dc
		// 
		// The input and output types are double.

		return tgamma(n + 1);
	}

	size_t nChoosek(const size_t n, size_t k) {
		// Computes the binomial coefficient (n choose k), i.e. the
		//  number of ways to choose k objects from a set of n objects.
		//
		// This code has been taken directly from https://stackoverflow.com/a/9331125
		//
		// See the following link for an explanation: https://stackoverflow.com/a/42285958

		if (k > n) return 0;
		if (k * 2 > n) k = n - k;
		if (k == 0) return 1;

		unsigned int result = n;
		for (unsigned int i = 2; i <= k; ++i) {
			result *= (n - i + 1);
			result /= i;
		}
		return result;
	}

	double nChoosek_gamma(const double n, const double k) {
		// Computes the binomial coefficient (n choose k), i.e. the
		// number of ways to choose k objects from a set of n objects.
		// 
		// This is a continuous version of nChoosek, extended by the use
		// of the gamma function.
		//
		// See: https://en.wikipedia.org/wiki/Binomial_coefficient#In_programming_languages

		if (k < 0) return 0;

		return exp(lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1));
	}

	std::vector<double> integrate_func_SISO(const std::vector<double>& ip_list,
		SISO_scalar_function func) {
		// 1-D numerical integration (trapezium rule) of a Single-Input Single-Output (SISO)
		// function.
		//
		// This function returns a vector containing the values of the integral
		// at all points in ip_list minus the value of the integral at the initial point of
		// ip_list

		std::vector<double> result;
		double h;

		result.push_back(0);

		for (size_t i = 1; i < ip_list.size(); i++) {
			h = ip_list[i] - ip_list[i - 1];
			result.push_back(result[i - 1] + 0.5 * h *
				(func(ip_list[i]) + func(ip_list[i - 1])));
		}

		return result;
	}

	template <typename Dist, typename Poly>
	double integrate_product_dist_polys(const std::vector<double>& ip_list,
		std::vector<std::shared_ptr<Dist>> ip_distributions,
		std::vector<std::shared_ptr<Poly>> ip_orth_polys) {
		// 1-D numerical integration (trapezium rule) of a the PDFs of a set of distributions
		// and a set of orthogonal polynomials.
		//
		// This function returns a vector containing the values of the integral
		// at all points in ip_list minus the value of the integral at the initial point of
		// ip_list
		
		std::vector<double> result;
		double h = 1;
		double current_eval = 1;
		double prev_eval = 1;

		for (size_t i = 0; i < ip_distributions.size(); i++) {
			prev_eval *= ip_distributions[i]->evaluate_PDF(ip_list[0]);
		}

		for (size_t i = 0; i < ip_orth_polys.size(); i++) {
			prev_eval *= ip_orth_polys[i]->evaluate(ip_list[0]);
		}

		result.push_back(0);

		for (size_t j = 1; j < ip_list.size(); j++) {
			h = ip_list[j] - ip_list[j - 1];
			current_eval = 1;

			for (size_t i = 0; i < ip_distributions.size(); i++) {
				current_eval *= ip_distributions[i]->evaluate_PDF(ip_list[j]);
			}

			for (size_t i = 0; i < ip_orth_polys.size(); i++) {
				current_eval *= ip_orth_polys[i]->evaluate(ip_list[j]);
			}
			
			result.push_back(result[j - 1] + 0.5 * h *
				(current_eval + prev_eval));

			prev_eval = current_eval;
		}

		return result.back();
	}

	// Template function instantiations
	template
		int factorial(const int n);

	template
		size_t factorial(const size_t n);

	template
	double integrate_product_dist_polys(const std::vector<double>& ip_list,
		std::vector<std::shared_ptr<NormalCDistribution>> ip_distributions,
		std::vector<std::shared_ptr<HermitePoly>> ip_orth_polys);

	template
	double integrate_product_dist_polys(const std::vector<double>& ip_list,
		std::vector<std::shared_ptr<UniformCDistribution>> ip_distributions,
		std::vector<std::shared_ptr<LegendrePoly>> ip_orth_polys);

	template
	double integrate_product_dist_polys(const std::vector<double>& ip_list,
		std::vector<std::shared_ptr<GammaCDistribution>> ip_distributions,
		std::vector<std::shared_ptr<LaguerrePoly>> ip_orth_polys);

	template
		double integrate_product_dist_polys(const std::vector<double>& ip_list,
		std::vector<std::shared_ptr<BetaCDistribution>> ip_distributions,
		std::vector<std::shared_ptr<JacobiPoly>> ip_orth_polys);
}
