#include <catch.hpp>

#include <cmath>
#include "finite_difference.hpp"

#include <iostream>

/*
Test_finitedifference
 
Test case for the finite difference function
 
@author: Matt Marti
@date: 2019-09-04
*/

TEST_CASE("Make sure function works at different orders of accuracy", "[finite_difference]") {
    using namespace std;
    using namespace Numerics;

    // Function input
    double h = .1;
    int N = 101;
    vector<double> thist = vector<double>(N);
    vector<double> yhist = vector<double>(N);
    vector<double> ydottruthhist = vector<double>(N);
    for (int ii = 0; ii < N; ii++) {
        double t = h * ii;
        thist[ii] = t;
        yhist[ii] = -9.81 / 2 * pow(t, 2) + 50 * t + 200;
        ydottruthhist[ii] = -9.81*t + 50;
    }

    // Test Function call
    vector<double> errvec = { 5e-1, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 2e-11, 2.5e-11, 5e-11 };
    for (size_t n = 1; n < errvec.size(); n++) {
        vector<double> ydothist = finite_difference(&yhist, &thist, n);
        double maxerr = -1.0;
        for (size_t ii = 0; ii < ydothist.size(); ii++) {
            double err = abs(ydothist[ii] - ydottruthhist[ii]);
            if (err > maxerr) {
                maxerr = err;
            }
        }
        REQUIRE(maxerr < errvec[n-1]);
    }
}

/*
TEST_CASE("Make sure function works for multiple dimensions of data", "[finite_difference]") {
    using namespace std;
    using namespace Numerics;

    // Function input
    h = .1;
    thist = (0:h:10)';
    yhist = [sin(thist), cos(thist)];
    n = 5;

    // Truth value
    ydottruthhist = [cos(thist), -sin(thist)];

    // Function call
    [ydothist] = finitedifference(yhist, thist, n);

    // Test Function call
    vector<double> errvec = { 5e-2; 5e-3; 2.5e-4; 2e-5; 2e-6; 5e-7; 2e-8; 2e-9; 1e-10; 1e-11; 2e-12 };
    for (int n = 1; n < errvec.size(); n++) {
        ydothist = finitedifference(yhist, thist, n);
        errhist = ydothist - ydottruthhist;
        maxerr = max(max(abs(errhist)));
        REQUIRE(maxerr < errvec(n));
    }
}
*/