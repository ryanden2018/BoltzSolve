#include <cstdio>
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
    printf("%20.19e\n",const_sqrt(1000.0));
    return 0;
}