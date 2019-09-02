#include <Cubic_Spline.hpp>
//#include <vector>

Numerics::Cubic_Spline(
    std::vector<double> xkvec,
    std::vector<double> fkvec,
    std::vector<double> fslope) {

}

double Numerics::Cubic_Spline::operator()(double xinter) {
    return 0;
}

std::vector<double> Numerics::Cubic_Spline::operator()(std::vector<double> xinter) {
    return std::vector<double>(10);
}