#ifndef CONSTMATHHPP
#define CONSTMATHHPP

constexpr double CONSTMATH_PI = 3.141592653589793238462;

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
    while(x > CONSTMATH_PI) x -= 2.0*CONSTMATH_PI;
    while(x < -CONSTMATH_PI) x += 2.0*CONSTMATH_PI;
    if(-CONSTMATH_PI/2.0 < x && x < CONSTMATH_PI/2.0) return const_sin_helper(x);
    if(x > CONSTMATH_PI/2.0) return const_sin_helper(CONSTMATH_PI-x);
    return const_sin_helper(-CONSTMATH_PI-x);
}

constexpr double const_cos(double x)
{
    return const_sin(CONSTMATH_PI/2.0 - x);
}

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

constexpr double exp_of_one = const_exp_helper(0.5)*const_exp_helper(0.5);

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

constexpr double const_norm(double x, double y)
{
    return const_sqrt(x*x+y*y);
}

#endif
