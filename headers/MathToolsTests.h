#ifndef H_MATH_TOOLS_TESTS
#define H_MATH_TOOLS_TESTS

namespace TestSuite {

	bool integrator_test();

	bool normal_distribution_test();

	bool unif_distribution_test_parametric();

	bool unif_distribution_test_mean();

	bool gamma_distribution_test_parametric();

	bool gamma_distribution_test_mean();

	bool beta_distribution_test_parametric();

	bool beta_distribution_test_mean();

	bool Legendre_Poly_Test();

	bool Hermite_Poly_Test();

	bool Laguerre_Poly_Test();

	bool Jacobi_Poly_Test();

	bool product_integrator_test();

	bool orthogonal_product_legendre_test();

	bool orthogonal_product_hermite_test();

	bool orthogonal_product_laguerre_test();

	bool galerkin_projection_test();
}

#endif
