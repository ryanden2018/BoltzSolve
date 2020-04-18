#ifndef HERMITE_HPP
#define HERMITE_HPP

#include <cmath>

const double HERMITE_PI = 3.141592653589793238462;

const double HERMITE_PI_QUARTER = std::pow(HERMITE_PI,0.25);

double hermiteEval(double x, int n)
{
    return n <= 0 
        ?
        std::exp(-0.5*x*x) / HERMITE_PI_QUARTER
        :
        (
            n == 1
            ?
                std::sqrt(2.0) * x * std::exp(-0.5*x*x) / HERMITE_PI_QUARTER
            :
                std::sqrt(2.0/n) * x * hermiteEval(x,n-1) - std::sqrt((n-1.0)/(1.0*n)) * hermiteEval(x,n-2)
        );
}

double hermiteDerivEval(double x, int n)
{
    return n <= 0
    ?
    (-1.0) * x * std::exp(-0.5 * x * x) * / HERMITE_PI_QUARTER
    :
    std::sqrt(n/2.0)*hermiteEval(x,n-1)-std::sqrt((n+1.0)/2.0)*hermiteEval(x,n+1);
}

#endif
