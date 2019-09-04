#ifndef NUMERICS_FINITE_DIFFERENCE_HPP
#define NUMERICS_FINITE_DIFFERENCE_HPP

#include <vector>

namespace Numerics {

    /*
    Finite Difference derivative function

    Foward, Central, Backwards finite difference calculation of derivative.
    This function uses the Central Finite Difference method to compute the
    derivative for the given data vector.At the ends of the array, central
    difference doesn't work, so the forward difference method and backwards
    difference methods are used instead.
        
    Note that this function decides for you to use forward and bakwards
    differencing functions at either end of the dataset.This cannot be
    turned off or changed.
        
    @arg
    yhist - N x M double matrix
            Function value time history, where N is the length of the
            dataset and M is the number different things to take the
            derivative of.M is usually 1.
    h - double
              Time step
    n - double(optional)
        Order of finite difference
        
    @return
    ydothist - N x M double matrix
               Finite difference derivative time history

    @author: Matt Marti
    @date: 2019-09-04
    */
    std::vector<double> finite_difference(std::vector<double> * yhist, std::vector<double> * thist, int n = 3);

} // namespace Numerics

#endif // NUMERICS_FINITE_DIFFERENCE_HPP