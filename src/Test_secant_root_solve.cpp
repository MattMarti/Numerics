#include <catch.hpp>

#include <secant_root_solve.hpp>
#include <cmath>

double test_function(double x) {
    using namespace std;
    return sin(x);
}

TEST_CASE("Test secant_root_solve on a function", "[secant-root-solve]") {
    double a, b, truth, x;
    a = -1;
    b = 3;
    std::function<double(double)> func = test_function;
    truth = 0.0;
    x = numerics::secant_root_solve(func, a, b);
    REQUIRE( std::abs(x - truth) < 1e-12 );
}