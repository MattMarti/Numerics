#ifndef NUMERICS_SECANT_ROOT_SOLVE_HPP
#define NUMERICS_SECANT_ROOT_SOLVE_HPP

#include <functional>

namespace Numerics
{
    /*
    * Solves for the roots of a function using the secant method
    * @arg
    * f         - function
    *             Function handle to solve root for
    * a         - double
    *             Upper bound
    * b         - double
    *             Lower bound
    * max_iters - int (optional)
    *             Maximum number of iterations
    * epsilon   - double (optional)
    *             Minimum error stopping criteria (in difference)
    * @return
    * x         - double
    *             Function root
    * @author: Matt Marti, Wade Foster
    * @date: 2019-08-31
    */
    double secant_root_solve(
        std::function<double(double)> function,
        double lower_bound, double upper_bound,
        size_t max_iters = 1000, double epsilon = 1e-12);

} // namespace Numerics

#endif // NUMERICS_SECANT_ROOT_SOLVE_H