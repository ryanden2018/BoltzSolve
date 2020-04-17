#ifndef BOLTZQ_HPP
#define BOLTZQ_HPP

#include "Hermite.hpp"

constexpr double const_sin_helper(double x)
{
    double uu = 1.0;
    double ans = 0.0;
    for(int i = 1; i < 25; i++)
    {
        uu *= x / i;
        if(i%2 == 1)
        {
            ans += uu;
            uu *= -1.0;
        }
    }
    return ans;
}

constexpr double const_sin(double x)
{
    while(x > HERMITE_PI) x -= 2.0*HERMITE_PI;
    while(x < -HERMITE_PI) x += 2.0*HERMITE_PI;
    if(-HERMITE_PI/2.0 < x && x < HERMITE_PI/2.0) return const_sin_helper(x);
    if(x > HERMITE_PI/2.0) return const_sin_helper(HERMITE_PI-x);
    return const_sin_helper(-HERMITE_PI-x);
}

constexpr double const_cos(double x)
{
    return const_sin(HERMITE_PI/2.0 - x);
}

#endif 
