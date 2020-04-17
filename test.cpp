#include <cstdio>
constexpr double HERMITE_PI = 3.141592653589793238462;
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
constexpr double abs(double x)
{
    return x > 0 ? x : -x;
}

constexpr double const_sqrt(double x)
{
    double ans = 1.0;
    while(abs(ans*ans-x)>1.0e-14)
    {
        ans = ans - (ans*ans-x)/(2*ans);
    }
    return ans;
}
int main()
{
    printf("%20.19e\n",const_sin(19.0*HERMITE_PI + 0.15));
    return 0;
}