#include <math.h>
#include <fstream>
#include "../headers/MathTools.h"

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
		result.push_back( result[i - 1] + 0.5 * h * 
			(func(ip_list[i]) + func(ip_list[i - 1])) );
	}

	return result;
}
