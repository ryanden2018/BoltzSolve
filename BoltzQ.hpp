#ifndef BOLTZQ_HPP
#define BOLTZQ_HPP

#include "Hermite.hpp"
#include <cmath>

const double BOLTZQ_PI = 3.141592653589793238462;

double computeBoltzQ(int ir, int jr, int kr, bool hardSphere)
{
    double h = 0.2;
    int Ntheta = 35;
    double val = 0.0;
    for(double r = 0.0; r < 10.0; r += h)
    {
        for(double vsx = -10.0; vsx < 10.0; vsx += h)
        {
            for(double vsy = -10.0; vsy < 10.0; vsy += h)
            {
                for(int i = 0; i < Ntheta; i++)
                {
                    double theta = (2.0 * BOLTZQ_PI * i)/Ntheta;
                    double tau = std::cos(theta)*(r-vsx) + std::sin(theta)*(-vsy);
                    double abs_tau = tau > 0.0 ? tau : -tau;
                    double vpx = r - std::cos(theta)*tau;
                    double vpy = -std::sin(theta)*tau;
                    double vspx = vsx + std::cos(theta)*tau;
                    double vspy = vsy + std::sin(theta)*tau;
                    double kernel = hardSphere ? abs_tau : 1.0;
                    double prefactor = r < 0.5*h ? 0.5 : 1.0;
                    val += prefactor * kernel * hermiteEval(std::sqrt(vpx*vpx+vpy*vpy),2*ir) 
                        * hermiteEval(std::sqrt(vspx*vspx+vspy*vspy),2*jr) 
                        * hermiteEval(r,2*kr) * 2.0 * std::sqrt(2.0) 
                        * h*h*h /Ntheta;
                    val -= prefactor * kernel * hermiteEval(r,2*ir) 
                        * hermiteEval(std::sqrt(vsx*vsx+vsy*vsy),2*jr) 
                        * hermiteEval(r,2*kr) * 2.0 * std::sqrt(2.0) 
                        * h*h*h /Ntheta;
                }
            }
        }
    }
    return val;
}



#endif