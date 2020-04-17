#ifndef HERMITE_HPP
#define HERMITE_HPP

const double HERMITE_PI = 3.141592653589793238462;

constexpr double const_exp_helper(double x)
{
    double uu = 1.0;
    double ans = uu;
    for(int i = 1; i < 16; i++)
    {
        uu *= x / i;
        ans += uu;
    }
    return ans;
}

const double exp_of_one = const_exp_helper(0.5)*const_exp_helper(0.5);

constexpr double const_exp(double x)
{
    return x < -0.51
    ?
    const_exp(x+1.0)/exp_of_one
    :
    (
        x > 0.51
        ?
        const_exp(x-1.0)*exp_of_one
        :
        const_exp_helper(x)
    );
}

constexpr double const_abs(double x)
{
    return x > 0 ? x : -x;
}

constexpr double const_sqrt(double x)
{
    double ans = 1.0;
    while(const_abs(ans*ans-x)>1.0e-14)
    {
        ans = ans - (ans*ans-x)/(2*ans);
    }
    return ans;
}

constexpr double hermiteEval(double x, int n)
{
    return n <= 0 
        ?
        const_exp(-0.5*x*x) / const_sqrt(const_sqrt(HERMITE_PI))
        :
        (
            n == 1
            ?
                const_sqrt(2.0) * x * const_exp(-0.5*x*x) / const_sqrt(const_sqrt(HERMITE_PI))
            :
                const_sqrt(2.0/n) * x * hermiteEval(x,n-1) - const_sqrt((n-1.0)/(1.0*n)) * hermiteEval(x,n-2)
        );
}

#endif
