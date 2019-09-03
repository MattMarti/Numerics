#ifndef NUMERICS_CUBIC_SPLINE_HPP
#define NUMERICS_CUBIC_SPLINE_HPP

#include <Eigen>
#include <vector>
#include <tuple>

namespace Numerics {

    /*
    Cubic Spline class allows for interpolation between data points

    @author: Matt Marti
    @date: 2019-09-03
     */
    class Cubic_Spline {

        Eigen::MatrixXd _akmat, _bkmat, _ckmat, _dkmat;
        std::vector<double> _xkvec;
        bool _extrapolation_enabled_flag = 0;

    public:

        /*
        Cubic spline constructor
        Solves for the parameters that describe a cubic spline. This
        function can solve the parameters of a cubic spline given
        data from a vectorized system. For example, the three space
        position of a planet in orbit can be interpolated.

        @arg
        xkvec         - 1 x n double matrix
                        Independent variable data points
        fkvec         - m x n double matrix
                        Dependent variable data points
        fslope        - m x 2 double matrix (optional)
                        Function slope at boundary points

        @return
        splineDataMat - m x n x 5 double matrix
                         Spline coefficient data matrix. Organized by
                         input data dimension, known value points, and
                         coefficient.

        @author: Matt Marti
        @date: 2019-09-03
        */
        Cubic_Spline(
            const std::vector<double> * xkvec,
            const std::vector<double> * fkvec,
            const std::vector<double> * fslope);

        /*
        Cubic spline interpolation function
        Interpolates function values at specified points using data
        from a solved cubic spline.

        @arg
        splineDataMat - m x n x 5 double matrix
                        Spline coefficient data matrix. Organized by
                        input data dimension, known value points, and
                        coefficient.
        xinter        - 1 x n double matrix
                        Interpolation points
        dflag         - bool (optional)
                        Optional flag to make the function return the
                        derivative.
                        False by default.

        @return
        finter        - n x 1 double matrix
                        Interpolated function value
        dfinter       - n x 1 double matrix
                        Interpolated function derivative value

        @author: Matt Marti
        @date: 2019-09-03
        */
        double operator()(double xinter);
        std::vector<double> operator()(std::vector<double>* xinter);
    };

}; // namespace Numerics

#endif // NUMERICS_CUBIC_SPLINE_HPP