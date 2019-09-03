#include <catch.hpp>
#include <Cubic_Spline.hpp>

/*
Test_cubicspline.m
 
Test case for the cubic spline function. Based on the driver script for
the solution to AOE 4404 Assignment 5, problems 1 to 3. Cubic spline and
Lagrange interpolation.

@author: Matt Marti
@date: 2019-09-01
*/

TEST_CASE("Test Cubic_Spline on a 1-D case", "[Cubic-Spline]") {
    using namespace std;

    // Given
    //f = @(x) sin(x);
    //df = @(x) cos(x);
    vector<double> xkvec, fkvec, xinterp;
    vector<vector<double>> fslope;
    for (int ii = 0; ii <= 20; ii++) xkvec[ii] = 
    xkvec = linspace(0, 10, 20);
    fkvec = f(xkvec);
    xinter = linspace(0, 10, 1000);
    fslope = [cos(xkvec(1)), cos(xkvec(end))]; % Clambed B.C.s

    // Run function
    cs = cubicspline(xkvec, fkvec, fslope);
    [finter, dfinter] = cs.interp(xinter, true);


    // Test Function truth values
    fitrue = f(xinter);
    error = fitrue - finter;
    maxerr = max(abs(error));
    assert(maxerr < 2.5e-4, 'Spline error too high');

    // Test Derivative truth values
    dfitrue = df(xinter);
    errord = dfitrue - dfinter;
    maxerrd = max(abs(errord));
    assert(maxerrd < 1.5e-3, 'Spline error too high');
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
    fitrue = f(xinter);
    error = fitrue - finter;
    maxerr = max(max(abs(error)));
    assert(maxerr < 2.5e-4, 'Spline error too high');

    // Test Derivative truth values
    dfitrue = df(xinter);
    errord = dfitrue - dfinter;
    maxerrd = max(max(abs(errord)));
    assert(maxerrd < 1.5e-3, 'Spline error too high');
}
*/