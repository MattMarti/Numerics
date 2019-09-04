#include <../eigen-git-mirror/Eigen/Eigen>
#include <cmath>
#include <finite_difference.hpp>

std::vector<double> Numerics::finite_difference(std::vector<double> * yhist, std::vector<double> * thist, int n) {
    using namespace std;
    using namespace Eigen;

    // Input checking
    int N = yhist->size();
    assert(N >= n + 1, "Not enough data points for given order");

    // Preallocate
    vector<double> ydothist = vector<double>(yhist->size());
    MatrixXd H = MatrixXd::Zero(n, n);
    VectorXd diffvec = VectorXd::Zero(n);
    VectorXd ydoti = VectorXd::Zero(n);
    double hj;
    int d, ii, jj, kk;

    // Forward difference method
    for (ii = 0; ii < n; ii++) {

        // Delta t matrix
        for (jj = 1; jj < n + 1; jj++) {
            d = 1;
            hj = (*thist)[ii + jj] - (*thist)[ii];
            diffvec(jj-1) = (*yhist)[ii + jj] - (*yhist)[ii];
            for (kk = 1; kk < n + 1; kk++) {
                d = d * kk;
                H(jj-1, kk-1) = pow(hj, kk) / d;
            }
        }

        // Compute derivative
        ydoti = H.partialPivLu().solve(diffvec);

        // Assign output
        ydothist[ii] = ydoti(0);
    }

    // Central difference method
    for (ii = n; ii < N - n; ii++) {

        // Delta t matrix
        for (jj = 1; jj < n + 1; jj++) {
            d = 1;
            hj = (*thist)[ii + jj] - (*thist)[ii];
            diffvec(jj-1) = (*yhist)[ii + jj] - (*yhist)[ii];
            for (kk = 1; kk < n + 1; kk++) {
                d = d * kk;
                H(jj-1, kk-1) = pow(hj, kk) / d;
            }
        }

        // Forward finite difference
        ydoti = 0.5 * H.partialPivLu().solve(diffvec);

        // Delta t matrix
        for (jj = 1; jj < n + 1; jj++) {
            d = 1;
            hj = (*thist)[ii] - (*thist)[ii - jj];
            diffvec(jj-1) = (*yhist)[ii] - (*yhist)[ii - jj];
            for (kk = 1; kk < n + 1; kk++) {
                d = d * kk;
                H(jj-1, kk-1) = pow(hj, kk) / d;
            }
        }

        // Compute derivative
        ydoti += 0.5 * H.partialPivLu().solve(diffvec);

        // Assign output
        ydothist[ii] = ydoti(0);
    }

    // Backwards difference method
    for (ii = N - n; ii < N; ii++) {

        // Delta t matrix
        for (jj = 1; jj < n + 1; jj++) {
            d = 1;
            hj = (*thist)[ii] - (*thist)[ii - jj];
            diffvec(jj-1) = (*yhist)[ii] - (*yhist)[ii - jj];
            for (kk = 1; kk < n + 1; kk++) {
                d = d * kk;
                H(jj - 1, kk - 1) = pow(hj, kk) / d;
            }
        }

        // Compute derivative
        ydoti = H.partialPivLu().solve(diffvec);

        // Assign output
        ydothist[ii] = ydoti[0];
    }
    
    return ydothist;
}