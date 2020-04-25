#include "Hermite.hpp"
#include <cstdio>

constexpr double INTEGRATE_PI = 3.141592653589793238462;

int main()
{
    for(int k = 0; k < 11; k++)
    {
        double val = 0.0;
        double h = 0.00001;
        for(double r = -10.0; r < 10.0; r += h)
        {
            double absr = r > 0.0 ? r : -r;
            val += 2.0 * INTEGRATE_PI * hermiteEval(r,2*k) * absr * h;
        }
        std::printf("%20.19e\n",val);
    }
    std::printf("------------\n");
    for(int k = 0; k < 11; k++)
    {
        double val = 0.0;
        double h = 0.00001;
        for(double r = -10.0; r < 10.0; r += h)
        {
            double absr = r > 0.0 ? r : -r;
            val += 2.0 * INTEGRATE_PI * hermiteEval(r,2*k) * absr * absr * absr * h;
        }
        std::printf("%20.19e\n",val);
    }
    return 0;
}
