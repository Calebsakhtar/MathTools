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

	// Given two doubles, this function returns whether they are within tol of each other.
	//
	// The input tol is optional, and has a default value of 1e-6.
	bool is_same_double(const double A, const double B, const double tol = 1e-6);

	// Prints a list of scalars into a ".csv" file
	void print_scalar_list(const std::vector<double>& scalar_list,
		const std::string& filename);

	// Prints a list of vectors into a ".csv" file
	void print_vector_list(const std::vector<Eigen::VectorXd>& vect_list,
		const std::string& filename);
	
	// Returns the factorial of an integer-like input.
	// The type of the output matches that of the input.
	template <typename int_type>
	int_type factorial(const int_type n);

	// Returns the value of the gamma function at the input plus one.
	// For positive integers, this is equivalent to computing the factorial.
	//
	// See for more details: https://www.cantorsparadise.com/the-beautiful-gamma-function-and-the-genius-who-discovered-it-8778437565dc
	// 
	// The input and output types are double.
	double factorial_gamma(const double n);

	// Computes the binomial coefficient (n choose k), i.e. the
	// number of ways to choose k objects from a set of n objects.
	//
	// This code has been taken directly from https://stackoverflow.com/a/9331125
	//
	// See the following link for an explanation: https://stackoverflow.com/a/42285958
	size_t nChoosek(const size_t n, size_t k);

	// Computes the binomial coefficient (n choose k), i.e. the
	// number of ways to choose k objects from a set of n objects.
	// 
	// This is a continuous version of nChoosek, extended by the use
	// of the gamma function.
	//
	// See: https://en.wikipedia.org/wiki/Binomial_coefficient#In_programming_languages
	double nChoosek_gamma(const double n, const double k);

	// 1-D numerical integration (trapezium rule) of a Single-Input Single-Output (SISO)
	// function.
	//
	// This function returns a vector containing the values of the integral
	// at all points in ip_list minus the value of the integral at the initial point of
	// ip_list
	std::vector<double> integrate_func_SISO(const std::vector<double>& ip_list,
		SISO_scalar_function func);

	// 1-D numerical integration (trapezium rule) of a the PDFs of a set of distributions
	// and a set of orthogonal polynomials.
	//
	// This function returns the total value of the integral
	template <typename Dist, typename Poly>
	double integrate_product_dist_polys(const std::vector<double>& ip_list,
		std::vector<std::shared_ptr<Dist>> ip_distributions,
		std::vector<std::shared_ptr<Poly>> ip_orth_polys);

	//// N-D numerical integration of a Multiple-Input Single-Output (MISO) function.
	//double integrate_func_MISO(const std::vector<Eigen::VectorXd>& ip_list,
	//	MISO_scalar_function func);
}

#endif
