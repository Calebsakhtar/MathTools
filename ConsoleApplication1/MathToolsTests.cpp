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

int main(){
    //integrator_test();
    //normal_distribution_test();
    //unif_distribution_test_parametric();
    //unif_distribution_test_mean();
    //gamma_distribution_test_parametric();
    //gamma_distribution_test_mean();
    beta_distribution_test_parametric();
    //beta_distribution_test_mean();
}
