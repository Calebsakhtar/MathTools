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

	size_t nChoosek(size_t n, size_t k)	{
		// Computes the binomial coefficient (n choose k), i.e. the
		//  number of ways to choose k objects from a set of n objects.
		//
		// This code has been taken directly from https://stackoverflow.com/a/9331125
		//
		// See the following link for an explanation: https://stackoverflow.com/a/42285958

		if (k > n) return 0;
		if (k * 2 > n) k = n - k;
		if (k == 0) return 1;

		int result = n;
		for (int i = 2; i <= k; ++i) {
			result *= (n - i + 1);
			result /= i;
		}
		return result;
	}
}
