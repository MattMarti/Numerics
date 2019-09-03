#ifndef CATCH_CONFIG_MAIN
#define CATCH_CONFIG_MAIN
#endif
#include <catch.hpp>

#include <secant_root_solve.hpp>
#include <cmath>

// test function
namespace Numerics {
    double _secant_root_solve_test_function(double x) {
        return std::sin(x);
    }
}

TEST_CASE("Test secant_root_solve on a function", "[secant-root-solve]") {
    double a, b, truth, x;
    a = -1;
    b = 3;
    truth = 0.0;
    std::function<double(double)> func = Numerics::_secant_root_solve_test_function;
    x = Numerics::secant_root_solve(func, a, b);
    REQUIRE( std::abs(x - truth) < 1e-12 );
}