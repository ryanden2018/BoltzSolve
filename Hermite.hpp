#ifndef HERMITE_HPP
#define HERMITE_HPP
#include <cmath>
const double HERMITE_PI = 3.141592653589793238462;

const double HERMITE_PI_QUARTER = std::sqrt(std::sqrt(HERMITE_PI));

double hermiteEval(double x, int n)
{
    double a = std::exp(-0.5*x*x) / HERMITE_PI_QUARTER;
    double b = std::sqrt(2.0) * x * std::exp(-0.5*x*x) / HERMITE_PI_QUARTER;
    for(int i = 0; i < n; i++)
    {
        double c =  std::sqrt(2.0) * std::pow((i+2.0),-0.5) * x * b - std::pow(i+1.0,0.5) * std::pow(i+2.0,-0.5) * a;
        a = b;
        b = c;
    }
    return a;
}

double hermiteDerivEval(double x, int n)
{
    return n <= 0
    ?
    (-1.0) * x * std::exp(-0.5 * x * x) / HERMITE_PI_QUARTER
    :
    std::sqrt(n/2.0)*hermiteEval(x,n-1)-std::sqrt((n+1.0)/2.0)*hermiteEval(x,n+1);
}


#endif