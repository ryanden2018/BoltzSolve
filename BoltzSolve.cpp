#include <cmath>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <SDL/SDL.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>

const double PI = 3.141592653589793238462;
const double N = 100000;


struct Velocity
{
    double vx;
    double vy;
    double vz;
    Velocity(double vx,double vy,double vz) : vx(vx),vy(vy),vz(vz) {}
};

struct ImpactParameter
{
    double omega_x;
    double omega_y;
    double omega_z;
    ImpactParameter(double omega_x,double omega_y,double omega_z)
        : omega_x(omega_x),omega_y(omega_y),omega_z(omega_z) {}
};

ImpactParameter RandomImpactParameter()
{
    double omega_x = 0.0;
    double omega_y = 0.0;
    double omega_z = 0.0;
    while(true)
    {
        double a = 2.0*((double)std::rand())/RAND_MAX - 1.0;
        double b = 2.0*((double)std::rand())/RAND_MAX - 1.0;
        if(a*a+b*b > 1.0) continue;
        double u = std::sqrt(1.0-a*a-b*b);
        omega_x = 2.0*a*u;
        omega_y = 2.0*b*u;
        omega_z = 1.0-2.0*(a*a+b*b);
        break;
    }
    ImpactParameter omega(omega_x,omega_y,omega_z);
    return omega;
}

void Evolve(std::vector<Velocity>& velocities)
{
    while(true)
    {
        int i = std::rand() % velocities.size();
        int j = std::rand() % velocities.size();
        if(i == j) continue;
        ImpactParameter omega = RandomImpactParameter();
        double vx = velocities[i].vx;
        double vy = velocities[i].vy;
        double vz = velocities[i].vz;
        double ux = velocities[j].vx;
        double uy = velocities[j].vy;
        double uz = velocities[j].vz;
        double tau = omega.omega_x*(vx-ux)+omega.omega_y*(vy-uy)+omega.omega_z*(vz-uz);
        velocities[i].vx = vx - omega.omega_x*tau;
        velocities[i].vy = vy - omega.omega_y*tau;
        velocities[i].vz = vz - omega.omega_z*tau;
        velocities[j].vx = ux + omega.omega_x*tau;
        velocities[j].vy = uy + omega.omega_y*tau;
        velocities[j].vz = uz + omega.omega_z*tau;
        break;
    }
}


extern "C"
{

std::vector<Velocity> velocities;

SDL_Surface *screen;

const int NUMPIXELS = 693;

void repaint()
{
	if (SDL_MUSTLOCK(screen)) SDL_LockSurface(screen);
	for (int i = 0; i < NUMPIXELS; i++) {
		for (int j = 0; j < NUMPIXELS; j++) {
			*((Uint32*)screen->pixels + i * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 255, 255, 255, 255);
		}
	}
    std::vector<double> vals;
    int K = 60;
    int num = velocities.size();
    for(int j = 0; j < K; j++)
    {
        double val = 0.0;
        double vx0 = 3.0*(2.0*j * (1.0/K) - 1.0);
        double vx1 = vx0 + 3.0*2.0/K;
        for(int i = 0; i < num; i++)
        {
            if(vx0 <= velocities[i].vx && velocities[i].vx<vx1) val += 1.0/num;
        }
        vals.push_back(val);
    }
    for (int j = 0; j < NUMPIXELS; j++)
    {   
        int idx = (j*(K-1))/NUMPIXELS ;
        double val = 10*vals[idx];
        int i = (int) (((val+1.0)/5.0)*NUMPIXELS);
        int ii = (int) (((1.0)/5.0)*NUMPIXELS);
        i = NUMPIXELS - i;
        ii = NUMPIXELS - ii;
        if(i < 0) i = 0;
        if(i > NUMPIXELS) i = NUMPIXELS;
        *((Uint32*)screen->pixels + i * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 0, 0, 0, 255);
        *((Uint32*)screen->pixels + ii * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 0, 0, 0, 255);
    }
	if (SDL_MUSTLOCK(screen)) SDL_UnlockSurface(screen);
	SDL_Flip(screen); 
}

int n = 0;
void sdl_init()
{
	if(n < 5)
	{
        n++;
		return;
	}
    repaint();

    for(int i = 0; i < 5000; i++) Evolve(velocities);
}

int main(int argc, char ** argv)
{
    std::srand(std::time(NULL));
	SDL_Init(SDL_INIT_VIDEO);
	screen = SDL_SetVideoMode(NUMPIXELS, NUMPIXELS, 32, SDL_SWSURFACE);

    for(int i = 0; i < N; i++)
    {
        ImpactParameter im = RandomImpactParameter();
        Velocity v(im.omega_x,im.omega_y,im.omega_z);
        velocities.push_back(v);
    }

	// emscripten_set_click_callback("canvas", 0, 1, mouseclick_callback);
	// emscripten_set_mousedown_callback("canvas", 0, 1, mousedown_callback);
	// emscripten_set_mouseup_callback("canvas", 0, 1, mouseup_callback);
	// emscripten_set_mousemove_callback("canvas", 0, 1, mousemove_callback);
	// emscripten_set_mouseleave_callback("canvas", 0, 1, mouseleave_callback);
	// emscripten_set_touchstart_callback("canvas", 0, 1, touchstart_callback);
	// emscripten_set_touchend_callback("canvas", 0, 1, touchend_callback);
	// emscripten_set_touchcancel_callback("canvas", 0, 1, touchcancel_callback);
	// emscripten_set_touchmove_callback("canvas", 0, 1, touchmove_callback);


	emscripten_set_main_loop(sdl_init,10,0);
}

}