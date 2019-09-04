#include <Eigen>
#include <cmath>
#include <Cubic_Spline.hpp>

Numerics::Cubic_Spline::Cubic_Spline(
        const std::vector<double> * xkvec,
        const std::vector<double> * fkvec,
        const std::vector<double> * fslope,
        bool extrapolation_enabled) {
    using namespace Eigen;
    using namespace std;

    // Extrapolation condition
    _extrapolation_enabled_flag = extrapolation_enabled;

    // Size
    size_t nx = xkvec->size();
    size_t m = 1;// fkvec.size();

    // Check input array sizes
    //size_t fkvec_size = fkvec->size();
    //size_t fslope_size = fslope->size();
    //static_assert(fkvec_size == nx, "Argument 'fkvec' does not match length of 'xkvec'");
    //static_assert(fslope_size == 2, "Argument 'fslope' is of incorrect size");
    bool boundary_is_clamped = 1;
    
    // Preallocate data storage
    _akmat = MatrixXd::Zero(nx, m);
    _bkmat = MatrixXd::Zero(nx, m);
    _ckmat = MatrixXd::Zero(nx, m);
    _dkmat = MatrixXd::Zero(nx, m);
    VectorXd xstar = VectorXd::Zero(nx);
    _xkvec = vector<double>(nx);
    for (size_t ii = 0; ii < nx; ii++) _xkvec[ii] = (*xkvec)[ii];

    // Build tri - diagonal system of equations
    size_t n = nx - 1;
    MatrixXd H = MatrixXd::Zero(nx, nx);
    VectorXd hkvec = VectorXd::Zero(nx);
    
    for (size_t ii = 0; ii < n; ii++) hkvec(ii) = _xkvec[ii + 1] - _xkvec[ii];
    H(0, 0) = 1;
    H(n, n) = 1;
    for (size_t k = 1; k < n; k++ ) {
        H(k, k - 1) = hkvec(k - 1);
        H(k, k) = 2 * (hkvec(k - 1) + hkvec(k));
        H(k, k + 1) = hkvec(k);
    }
    if (boundary_is_clamped) {
        H(0, 0) = 2 * hkvec(0);
        H(0, 1) = hkvec(0);
        H(n, n-1) = hkvec(n-1);
        H(n, n) = 2 * hkvec(n-1);
    }

    // Generate right hand side
    for (size_t ii = 0; ii < nx; ii++) _akmat(ii, 0) = (*fkvec)[ii];
    for (size_t k = 1; k < n; k++) {
        xstar(k) = 3 * ((_akmat(k + 1, 0) - _akmat(k, 0)) / hkvec(k) \
            - (_akmat(k, 0) - _akmat(k - 1, 0)) / hkvec(k - 1));
    }
    if (boundary_is_clamped) {
        xstar(0) = 3 * ( (_akmat(1, 0) - _akmat(0, 0)) / hkvec(0) - (*fslope)[0] );
        xstar(n) = 3 * ( (*fslope)[1] - (_akmat(n, 0) - _akmat(n-1, 0)) / hkvec(n-1) );
    }

    // Solve tri - diagonal system of equations
    _ckmat = H.partialPivLu().solve(xstar);

    // Compute bkmat and dkmat
    for (size_t k = 0; k < n; k++) {
        _bkmat(k, 0) = (_akmat(k + 1, 0) - _akmat(k, 0)) / hkvec(k) \
                     - hkvec(k) * (2 * _ckmat(k, 0) + _ckmat(k + 1, 0)) / 3;
        _dkmat(k, 0) = (_ckmat(k + 1, 0) - _ckmat(k, 0)) / (3 * hkvec(k));
    }
    if (boundary_is_clamped) {
        _bkmat(n, 0) = (*fslope)[1];
    }
}

double Numerics::Cubic_Spline::operator()(double xinter) {
    
    // Size
    size_t nx = _akmat.rows();
    size_t n = nx - 1;

    // Check that interpolated value is within function range
    bool is_within_range = xinter >= _xkvec[0] && xinter <= _xkvec[n];
    if (!_extrapolation_enabled_flag) assert(is_within_range);

    // Find x value just below xinter(i)
    size_t k = 1;
    if (is_within_range) {
        while (k < nx) {
            if (xinter < _xkvec[k]) {
                k = k - 1;
                break;
            }
            k = k + 1;
        }
        if (k >= nx) { // Point is on upper boundary) {
            return _akmat(nx, 0);
        }
    }
    else if (xinter > _xkvec[n]) {
        k = n;
    }
    else {
        k = 0;
    }

    // Spline interpolation
    double hi = xinter - _xkvec[k];
    return _akmat(k, 0) + _bkmat(k, 0)*hi + _ckmat(k, 0)*pow(hi, 2) + _dkmat(k, 0)*pow(hi, 3);
}

std::vector<double> Numerics::Cubic_Spline::operator()(std::vector<double>* xinter) {
    std::vector<double> finter = std::vector<double>(xinter->size());
    for (size_t ii = 0; ii < xinter->size(); ii++) {
        finter[ii] = Cubic_Spline::operator()((*xinter)[ii]);
    }
    return finter;
}