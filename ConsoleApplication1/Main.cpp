#include "../headers/MathToolsTests.h"
#include "../headers/MathTools.h"

double f(double x) {
    return x;
}

int main() {
    //TestSuite::integrator_test();
    //TestSuite::normal_distribution_test();
    //TestSuite::unif_distribution_test_parametric();
    //TestSuite::unif_distribution_test_mean();
    //TestSuite::gamma_distribution_test_parametric();
    //TestSuite::gamma_distribution_test_mean();
    //TestSuite::wa_gamma_distribution_test_parametric();
    //TestSuite::wa_gamma_distribution_test_mean();
    //TestSuite::beta_distribution_test_parametric();
    //TestSuite::beta_distribution_test_mean();
    //TestSuite::wa_beta_distribution_test_parametric();
    //TestSuite::wa_beta_distribution_test_mean();
    //TestSuite::legendre_poly_test();
    //TestSuite::hermite_poly_test();
    //TestSuite::laguerre_poly_test();
    //TestSuite::jacobi_poly_test();
    //TestSuite::product_integrator_test();
    //TestSuite::orthogonal_product_legendre_test();
    //TestSuite::orthogonal_product_hermite_test();
    //TestSuite::orthogonal_product_laguerre_test();
    //TestSuite::cdf_icdf_normal_test();
    //TestSuite::cdf_icdf_uniform_test();
    // TestSuite::cdf_icdf_gamma_test();
    TestSuite::galerkin_projection_test();
}
