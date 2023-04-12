#include "../headers/MathTools.h"

double f(double x) {
    return x;
}

void integrator_test() {
    std::vector<double> x = MathTools::linspace(0, 10, 101);
    std::vector<double> x_integ_analytic = MathTools::linspace(0, 10, 101);
    std::vector<double> x_integ_numeric;

    for (size_t i = 0; i < x.size(); i++) {
        x_integ_analytic[i] = 0.5 * x_integ_analytic[i] * x_integ_analytic[i];
    }

    x_integ_numeric = MathTools::integrate_func_SISO(x, f);

    MathTools::print_scalar_list(x, "x");
    MathTools::print_scalar_list(x_integ_analytic, "xAnalytic");
    MathTools::print_scalar_list(x_integ_numeric, "xNumeric");
}

void normal_distribution_test() {
    std::vector<double> x = MathTools::linspace(2, 8, 201);
    std::vector<double> PDF_values;
    std::vector<double> samples;

    std::default_random_engine generator;
    MathTools::NormalCDistribution norm(5, 1);

    for (size_t i = 0; i < x.size(); i++) {
        PDF_values.push_back(norm.evaluate_PDF(x[i]));
        samples.push_back(norm.sample(generator));
    }

    MathTools::print_scalar_list(x, "xDistribution");
    MathTools::print_scalar_list(PDF_values, "PDF");
    MathTools::print_scalar_list(samples, "Samples");
}

void unif_distribution_test_parametric() {
    std::vector<double> x = MathTools::linspace(0, 6, 201);
    std::vector<double> PDF_values;
    std::vector<double> samples;

    std::default_random_engine generator;
    MathTools::UniformCDistribution unif(1, 5);

    for (size_t i = 0; i < x.size(); i++) {
        PDF_values.push_back(unif.evaluate_PDF(x[i]));
        samples.push_back(unif.sample(generator));
    }

    MathTools::print_scalar_list(x, "xDistribution");
    MathTools::print_scalar_list(PDF_values, "PDF");
    MathTools::print_scalar_list(samples, "Samples");
}

void unif_distribution_test_mean() {
    std::vector<double> x = MathTools::linspace(0, 6, 201);
    std::vector<double> PDF_values;
    std::vector<double> samples;

    std::default_random_engine generator;
    MathTools::UniformCDistribution unif;

    unif.set_params_mean_stdev(3, 4.0 / sqrt(12.0));

    for (size_t i = 0; i < x.size(); i++) {
        PDF_values.push_back(unif.evaluate_PDF(x[i]));
        samples.push_back(unif.sample(generator));
    }

    MathTools::print_scalar_list(x, "xDistribution");
    MathTools::print_scalar_list(PDF_values, "PDF");
    MathTools::print_scalar_list(samples, "Samples");
}

void gamma_distribution_test_parametric() {
    std::vector<double> x = MathTools::linspace(0, 6, 201);
    std::vector<double> PDF_values;
    std::vector<double> samples;

    std::default_random_engine generator;
    MathTools::GammaCDistribution gamma(1, 2);

    for (size_t i = 0; i < x.size(); i++) {
        PDF_values.push_back(gamma.evaluate_PDF(x[i]));
        samples.push_back(gamma.sample(generator));
    }

    MathTools::print_scalar_list(x, "xDistribution");
    MathTools::print_scalar_list(PDF_values, "PDF");
    MathTools::print_scalar_list(samples, "Samples");
}

void gamma_distribution_test_mean() {
    std::vector<double> x = MathTools::linspace(0, 6, 201);
    std::vector<double> PDF_values;
    std::vector<double> samples;

    std::default_random_engine generator;
    MathTools::GammaCDistribution gamma;

    gamma.set_params_mean_stdev(0.5, 0.5);

    for (size_t i = 0; i < x.size(); i++) {
        PDF_values.push_back(gamma.evaluate_PDF(x[i]));
        samples.push_back(gamma.sample(generator));
    }

    MathTools::print_scalar_list(x, "xDistribution");
    MathTools::print_scalar_list(PDF_values, "PDF");
    MathTools::print_scalar_list(samples, "Samples");
}

void beta_distribution_test_parametric() {
    std::vector<double> x = MathTools::linspace(0, 1, 201);
    std::vector<double> PDF_values;
    std::vector<double> samples;

    std::default_random_engine generator;
    MathTools::BetaCDistribution beta(2, 5);

    for (size_t i = 0; i < x.size(); i++) {
        PDF_values.push_back(beta.evaluate_PDF(x[i]));
        samples.push_back(beta.sample(generator));
    }

    MathTools::print_scalar_list(x, "xDistribution");
    MathTools::print_scalar_list(PDF_values, "PDF");
    MathTools::print_scalar_list(samples, "Samples");
}

void beta_distribution_test_mean() {
    std::vector<double> x = MathTools::linspace(0, 1, 201);
    std::vector<double> PDF_values;
    std::vector<double> samples;

    std::default_random_engine generator;
    MathTools::BetaCDistribution beta;

    beta.set_params_mean_stdev(0.2857142857142857, 0.15971914124998499);

    for (size_t i = 0; i < x.size(); i++) {
        PDF_values.push_back(beta.evaluate_PDF(x[i]));
        samples.push_back(beta.sample(generator));
    }

    MathTools::print_scalar_list(x, "xDistribution");
    MathTools::print_scalar_list(PDF_values, "PDF");
    MathTools::print_scalar_list(samples, "Samples");
}

void Legendre_Poly_Test() {
    std::vector<double> x = MathTools::linspace(-2, 2, 201);
    std::vector<double> P0, P1, P2, P3, P4, P5;

    MathTools::LegendrePoly Poly0(0), Poly1(1), Poly2(2), Poly3(3),
        Poly4(4), Poly5(5);
    
    for (size_t i = 0; i < x.size(); i++) {
        P0.push_back(Poly0.evaluate(x[i]));
        P1.push_back(Poly1.evaluate(x[i]));
        P2.push_back(Poly2.evaluate(x[i]));
        P3.push_back(Poly3.evaluate(x[i]));
        P4.push_back(Poly4.evaluate(x[i]));
        P5.push_back(Poly5.evaluate(x[i]));
    }

    MathTools::print_scalar_list(x, "xPolynomials");
    MathTools::print_scalar_list(P0, "Poly00");
    MathTools::print_scalar_list(P1, "Poly01");
    MathTools::print_scalar_list(P2, "Poly02");
    MathTools::print_scalar_list(P3, "Poly03");
    MathTools::print_scalar_list(P4, "Poly04");
    MathTools::print_scalar_list(P5, "Poly05");
}

void Hermite_Poly_Test() {
    std::vector<double> x = MathTools::linspace(-2, 2, 201);
    std::vector<double> P0, P1, P2, P3, P4, P5;

    MathTools::HermitePoly Poly0(0), Poly1(1), Poly2(2), Poly3(3),
        Poly4(4), Poly5(5);

    for (size_t i = 0; i < x.size(); i++) {
        P0.push_back(Poly0.evaluate(x[i]));
        P1.push_back(Poly1.evaluate(x[i]));
        P2.push_back(Poly2.evaluate(x[i]));
        P3.push_back(Poly3.evaluate(x[i]));
        P4.push_back(Poly4.evaluate(x[i]));
        P5.push_back(Poly5.evaluate(x[i]));
    }

    MathTools::print_scalar_list(x, "xPolynomials");
    MathTools::print_scalar_list(P0, "Poly00");
    MathTools::print_scalar_list(P1, "Poly01");
    MathTools::print_scalar_list(P2, "Poly02");
    MathTools::print_scalar_list(P3, "Poly03");
    MathTools::print_scalar_list(P4, "Poly04");
    MathTools::print_scalar_list(P5, "Poly05");
}

void Laguerre_Poly_Test() {
    std::vector<double> x = MathTools::linspace(-5, 15, 1001);
    std::vector<double> P0, P1, P2, P3, P4, P5;

    MathTools::LaguerrePoly Poly0(0), Poly1(1), Poly2(2), Poly3(3),
        Poly4(4), Poly5(5);

    for (size_t i = 0; i < x.size(); i++) {
        P0.push_back(Poly0.evaluate(x[i]));
        P1.push_back(Poly1.evaluate(x[i]));
        P2.push_back(Poly2.evaluate(x[i]));
        P3.push_back(Poly3.evaluate(x[i]));
        P4.push_back(Poly4.evaluate(x[i]));
        P5.push_back(Poly5.evaluate(x[i]));
    }

    MathTools::print_scalar_list(x, "xPolynomials");
    MathTools::print_scalar_list(P0, "Poly00");
    MathTools::print_scalar_list(P1, "Poly01");
    MathTools::print_scalar_list(P2, "Poly02");
    MathTools::print_scalar_list(P3, "Poly03");
    MathTools::print_scalar_list(P4, "Poly04");
    MathTools::print_scalar_list(P5, "Poly05");
}

void product_integrator_test() {
    // Values of the spectral variable eta
    const std::vector<double> ip_list = MathTools::linspace(-50, 50, 10001);

    // Other variables
    double result = 0; // Result value

    // Distribution shared pointers
    std::shared_ptr<MathTools::UniformCDistribution>
        germ_ptr(new MathTools::UniformCDistribution(0, 1)); // Germ Distribution

    // Orthogonal Polynomials
    std::shared_ptr<MathTools::LegendrePoly> current_poly_ptr(new MathTools::LegendrePoly(0));

    // Functions to be evaluated in the product integral
    std::vector<std::shared_ptr<MathTools::UniformCDistribution>>
        distributions = { germ_ptr };
    std::vector<std::shared_ptr<MathTools::LegendrePoly>> polys = { current_poly_ptr };

    // Evaluate the product integral
    result = MathTools::integrate_product_dist_polys(ip_list, distributions, polys);
}

void galerkin_projection_test() {
    // Values of the spectral variable eta
    const std::vector<double> ip_list = MathTools::linspace(-50, 50, 10001);

    // Other variables
    const size_t max_order = 2; // Maximum order of the orthogonal polynomials used
    double num = 0; // Numerator value
    double denom = 0; // Denominator value

    // Distribution shared pointers
    std::shared_ptr<MathTools::NormalCDistribution> 
        lambda_ptr(new MathTools::NormalCDistribution(7, 2)); // Uncertain Parameter Distribution
    std::shared_ptr<MathTools::NormalCDistribution> 
        germ_ptr(new MathTools::NormalCDistribution(0, 1)); // Germ Distribution

    // Orthogonal Polynomials
    std::shared_ptr<MathTools::HermitePoly> current_poly_ptr(new MathTools::HermitePoly(0));

    // Functions to be evaluated in the numerator of the Galerkin Projection 
    std::vector<std::shared_ptr<MathTools::NormalCDistribution>> 
        distributions_num = { lambda_ptr, germ_ptr };
    std::vector<std::shared_ptr<MathTools::HermitePoly>> polys_num = { current_poly_ptr };

    // Functions to be evaluated in the denominator of the Galerkin Projection
    std::vector<std::shared_ptr<MathTools::NormalCDistribution>>
        distributions_den = { germ_ptr };
    std::vector<std::shared_ptr<MathTools::HermitePoly>> polys_den = 
    { current_poly_ptr, current_poly_ptr };
    
    // Storage vectors for the coefficients of the uncertain variables
    std::vector<double> lambda_coeffs;
    std::vector<double> u_coeffs;

    // Loop to find the coefficients of Lambda
    for (size_t i = 0; i < max_order; i++) {
        // To save memory, only one polynomial is used, and is order is updated
        current_poly_ptr->set_order(i);

        // Evaluate the numerator and denominator of the Galerkin projection for Lambda
        num = MathTools::integrate_product_dist_polys(ip_list, distributions_num, polys_num);
        denom = MathTools::integrate_product_dist_polys(ip_list, distributions_den, polys_den);

        // Store the current coefficient of Lambda
        lambda_coeffs.push_back(num / denom);
    }
}

int main(){
    //integrator_test();
    //normal_distribution_test();
    //unif_distribution_test_parametric();
    //unif_distribution_test_mean();
    //gamma_distribution_test_parametric();
    //gamma_distribution_test_mean();
    //beta_distribution_test_parametric();
    //beta_distribution_test_mean();
    //Legendre_Poly_Test();
    //Hermite_Poly_Test();
    //Laguerre_Poly_Test();
    //product_integrator_test();
    //galerkin_projection_test();
}
