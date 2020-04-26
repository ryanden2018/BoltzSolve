#ifndef BOLTZQ_HPP
#define BOLTZQ_HPP

#include "Hermite.hpp"
#include <cmath>

const float BOLTZQ_PI = 3.141592653589793238462;

float computeBoltzQ(int ir, int jr, int kr, bool hardSphere)
{
    float h = 0.33*0.25;
    int Ntheta = 20*4;
    float val = 0.0;
    for(float r = 0.0; r < 10.0; r += h)
    {
        for(float vsx = -10.0; vsx < 10.0; vsx += h)
        {
            for(float vsy = -10.0; vsy < 10.0; vsy += h)
            {
                for(int i = 0; i < Ntheta; i++)
                {
                    float theta = (2.0 * BOLTZQ_PI * i)/Ntheta;
                    float tau = std::cos(theta)*(r-vsx) + std::sin(theta)*(-vsy);
                    float abs_tau = tau > 0.0 ? tau : -tau;
                    float vpx = r - std::cos(theta)*tau;
                    float vpy = -std::sin(theta)*tau;
                    float vspx = vsx + std::cos(theta)*tau;
                    float vspy = vsy + std::sin(theta)*tau;
                    float kernel = (hardSphere ? abs_tau : 1.0);
                    float prefactor = (r < 0.5*h || r > 10.0-h ? 0.5 : 1.0);
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
