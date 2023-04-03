#ifndef MATH_TOOLS_H
#define MATH_TOOLS_H

#include <vector>
#include <string>
#include"../eigen/Eigen/Dense"
#include "../headers/Distribution.h"
#include "../headers/OrthPoly.h"

typedef double(*SISO_scalar_function)(const double); // Single-Input Single-Output function
// typedef double(*MISO_scalar_function)(const Eigen::VectorXd&); // Multiple-Input Single-Output function

namespace MathTools {

	// Returns a std::vector containing num_points evenly spaced between start_val and end_val.
	//
	// num_points is inclusive of start_val and end_val
	std::vector<double> linspace(const double start_val, const double end_val,
		const size_t num_pts);

	// Prints a list of scalars into a ".csv" file
	void print_scalar_list(const std::vector<double>& scalar_list,
		const std::string& filename);

	// Prints a list of vectors into a ".csv" file
	void print_vector_list(const std::vector<Eigen::VectorXd>& vect_list,
		const std::string& filename);

	// 1-D numerical integration (trapezium rule) of a Single-Input Single-Output (SISO)
	// function.
	//
	// This function returns a vector containing the values of the integral
	// at all points in ip_list minus the value of the integral at the initial point of
	// ip_list
	std::vector<double> integrate_func_SISO(const std::vector<double>& ip_list,
		SISO_scalar_function func);

	// Computes the binomial coefficient (n choose k), i.e. the
	//  number of ways to choose k objects from a set of n objects.
	//
	// This code has been taken directly from https://stackoverflow.com/a/9331125
	//
	// See the following link for an explanation: https://stackoverflow.com/a/42285958
	size_t nChoosek(size_t n, size_t k);

	//// N-D numerical integration of a Multiple-Input Single-Output (MISO) function.
	//double integrate_func_MISO(const std::vector<Eigen::VectorXd>& ip_list,
	//	MISO_scalar_function func);
}

#endif
