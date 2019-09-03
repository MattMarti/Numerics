#include <catch.hpp>

#include <functional>
#include <Cubic_Spline.hpp>

/*
Test case for the cubic spline function. Based on the driver script for
the solution to AOE 4404 Assignment 5, problems 1 to 3. Cubic spline and
Lagrange interpolation.

@author: Matt Marti
@date: 2019-09-01
*/

TEST_CASE("Test Cubic_Spline on a 1-D case", "[Cubic-Spline]") {
    using namespace std;
    using namespace Numerics;

    // Given
    //f = @(x) sin(x);
    //df = @(x) cos(x);
    vector<double> xkvec, fkvec, xinter, ftrue;
    vector<double> fslope;
    xkvec = vector<double>(21);
    fkvec = vector<double>(21);
    for (int ii = 0; ii <= 21; ii++) {
        xkvec[ii] = 0.5*ii; // xkvec = linspace(0, 10, 20);
        fkvec[ii] = sin(xkvec[ii]); // fkvec = f(xkvec);
    }
    xinter = vector<double>(1001);
    ftrue = vector<double>(1001);
    for (int ii = 0; ii <= 1001; ii++) {
        xinter[ii] = 0.01*ii; // xinter = linspace(0, 10, 1000);
        ftrue[ii] = sin(xinter[ii]);
    }
    fslope = vector<double>(2); // fslope = [cos(xkvec(1)), cos(xkvec(end))]; % Clambed B.C.s
    fslope[0] = cos(xkvec[0]);
    fslope[1] = cos(xkvec[20]);

    // Run function
    Cubic_Spline cs = Cubic_Spline(&xkvec, &fkvec, &fslope);
    std::vector<double> finter = cs(&xinter);

    // Test Function truth values
    double err = INFINITY;
    double errii;
    for (size_t ii = 0; ii < finter.size(); ii++) {
        errii = abs(finter[ii] - ftrue[ii]);
        if (err > errii) err = errii;
    }
    REQUIRE(err < 2.5e-4);
}

TEST_CASE("Test for extrapolation", "[Cubic-Spline]") {

}

/*
TEST_CASE("Test Cubic_Spline on a Multi-Dimensional case", "[Cubic-Spline]") {
    using namespace std;

    // Given
    f = @(x)[sin(x); cos(x)];
    df = @(x)[cos(x); -sin(x)];
    xkvec = linspace(0, 10, 20);
    fkvec = f(xkvec);
    xinter = linspace(0, 10, 1000);
    fslope = [cos(xkvec(1)), cos(xkvec(end)); -sin(xkvec(1)), -sin(xkvec(end))]; % Clambed B.C.s

    // Run function
    cs = cubicspline(xkvec, fkvec, fslope);
    [finter, dfinter] = cs.interp(xinter, true);

    // Test Function truth values
    ftrue = f(xinter);
    error = ftrue - finter;
    maxerr = max(max(abs(error)));
    assert(maxerr < 2.5e-4, 'Spline error too high');

    // Test Derivative truth values
    dftrue = df(xinter);
    errord = dftrue - dfinter;
    maxerrd = max(max(abs(errord)));
    assert(maxerrd < 1.5e-3, 'Spline error too high');
}
*/