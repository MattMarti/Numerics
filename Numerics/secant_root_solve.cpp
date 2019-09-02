// secant_root_solve.cpp

#include <secant_root_solve.hpp>
#include <cmath>

double Numerics::secant_root_solve(
        std::function<double(double)> function,
        double lower_bound, double upper_bound,
        size_t max_iters, double epsilon) {

    // Check that root exists (Intermediate value theorem)
    double f_upper = function(upper_bound);
    double f_lower = function(lower_bound);
    // assert(f(upper_bound) * f(lower_bound) <= 0, 'Root does not exist');

    // Initialize loop
    double x_est = lower_bound;
    if (std::abs(f_upper) < std::abs(f_lower)) x_est = upper_bound;

    double error = std::abs(function(x_est));
    size_t iters = 0;
    for (; error > epsilon && iters < max_iters; ++iters)
    {
        double xi = upper_bound - f_upper * (lower_bound - upper_bound)
            / (f_lower - f_upper);

        // Assign this guess to next bounds
        double y = function(xi);
        if (y*f_upper < 0) lower_bound = xi;
        else if (y*f_lower < 0) upper_bound = xi;

        // Measure error
        error = std::abs((xi - x_est) / x_est);

        // Iterate
        x_est = xi;
    }

    return x_est;
}