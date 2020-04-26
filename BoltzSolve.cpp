#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <stdio.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <SDL/SDL.h>
#include <queue>
#include "Hermite.hpp"
#include <cstdio>
#include <cstdlib>
#include "BoltzQVals.hpp"
#include "integrals.hpp"

const int N = 11;

void computeOp(double *in, double* out)
{
    for(int k = 0; k < N; k++)
    {
        out[k] = 0.0;
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                out[k] += in[i] * boltzQMaxwell[i][j][k] * in[j];
            }
        }
    }
}

void scale(double *out, double lambda)
{
    for(int i = 0; i < N; i++)
    {
        out[i] *= lambda;
    }
}

void add(double *in, double *out)
{
    for(int i = 0; i < N; i++)
    {
        out[i] += in[i];
    }
}

void swap(double *out1, double *out2)
{
    double temp = 0.0;
    for(int i = 0; i < N; i++)
    {
        temp = out1[i];
        out1[i] = out2[i];
        out2[i] = temp;
    }
}

void init(double *out)
{
    for(int i = 0; i < N; i++)
    {
        out[i] = 0.0;
    }
}

double integral(double *in)
{
    double res = 0.0;
    for(int i = 0; i < N; i++)
    {
        res += hermiteIntegrals[i] * in[i];
    }
    return res;
}

double moment(double *in)
{
    double res = 0.0;
    for(int i = 0; i < N; i++)
    {
        res += hermiteMoments[i] * in[i];
    }
    return res;
}

const double PI = 3.141592653589793238462;

extern "C"
{

SDL_Surface *screen;

double *y;
double *k1;
double *k2;
double *k3;
double *k4;
double *tmp;
const int NUMPIXELS = 693;
double gamma = 0.0;
double lambda = 1.0;
double mass = 0.0;
double energy = 0.0;

void repaint()
{
	if (SDL_MUSTLOCK(screen)) SDL_LockSurface(screen);
	for (int i = 0; i < NUMPIXELS; i++) {
		for (int j = 0; j < NUMPIXELS; j++) {
			*((Uint32*)screen->pixels + i * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 255, 255, 255, 255);
		}
	}
    double num = 0.0;
    double num2 = 0.0;
    for (int j = 0; j < NUMPIXELS; j++)
    {
        double r = 20.0*(j-0.5*NUMPIXELS)/NUMPIXELS;
        double val = 0.0;
        for(int n = 0; n < N; n++)
        {
            val += y[n] * lambda * hermiteEval(lambda*r,2*n);
        }
        num += val*(r>0.0?r:-r);
        num2 += val*std::pow((r>0.0?r:-r),3);
        int i = (int) (((val+1.0)/5.0)*NUMPIXELS);
        int ii = (int) (((1.0)/5.0)*NUMPIXELS);
        i = NUMPIXELS - i;
        ii = NUMPIXELS - ii;
        if(i < 0) i = 0;
        if(i > NUMPIXELS) i = NUMPIXELS;
        *((Uint32*)screen->pixels + i * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 0, 0, 0, 255);
        *((Uint32*)screen->pixels + ii * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 0, 0, 0, 255);
    }
    printf("%7.6e\t%7.6e\n",integral(y)/mass,moment(y)/energy);
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



    double h = 0.00001;
    for(int i = 0; i < 100; i++)
    {
        computeOp(y,k1);
        init(tmp);
        add(k1,tmp);
        scale(tmp,h/2.0);
        add(y,tmp);
        computeOp(tmp,k2);
        init(tmp);
        add(k2,tmp);
        scale(tmp,h/2.0);
        add(y,tmp);
        computeOp(tmp,k3);
        init(tmp);
        add(k3,tmp);
        scale(tmp,h);
        add(y,tmp);
        computeOp(tmp,k4);
        init(tmp);
        add(k1,tmp);
        add(k2,tmp);
        add(k2,tmp);
        add(k3,tmp);
        add(k3,tmp);
        add(k4,tmp);
        scale(tmp,h/6.0);
        add(y,tmp);
        swap(y,tmp);
    }

	repaint();
}

int main(int argc, char ** argv)
{
    y = (double *)std::malloc(N*sizeof(double));
    k1 = (double *)std::malloc(N*sizeof(double));
    k2 = (double *)std::malloc(N*sizeof(double));
    k3 = (double *)std::malloc(N*sizeof(double));
    k4 = (double *)std::malloc(N*sizeof(double));
    tmp = (double *)std::malloc(N*sizeof(double));

    init(y);
    init(k1);
    init(k2);
    init(k3);
    init(k4);
    init(tmp);

    // y[0] = 1.0;
    // y[1] = 0.5;
    // y[2] = 0.3;
    // y[10] = 1.0;
    // mass = integral(y);
    // energy = moment(y);
    y[1] = 1.0;
    double m1 = integral(y);
    y[1] = 0.0;
    y[2] = 1.0;
    double m2 = integral(y);
    y[2] = 0.0;

    y[0] = 1.0;
    y[1] = 0.1*m2/m1;
    y[2] = -0.1;
    mass = integral(y);
    energy = moment(y);


    // std::free(y);
    // std::free(k1);
    // std::free(k2);
    // std::free(k3);
    // std::free(k4);
    // std::free(tmp);
	SDL_Init(SDL_INIT_VIDEO);
	screen = SDL_SetVideoMode(NUMPIXELS, NUMPIXELS, 32, SDL_SWSURFACE);

	// emscripten_set_click_callback("canvas", 0, 1, mouseclick_callback);
	// emscripten_set_mousedown_callback("canvas", 0, 1, mousedown_callback);
	// emscripten_set_mouseup_callback("canvas", 0, 1, mouseup_callback);
	// emscripten_set_mousemove_callback("canvas", 0, 1, mousemove_callback);
	// emscripten_set_mouseleave_callback("canvas", 0, 1, mouseleave_callback);
	// emscripten_set_touchstart_callback("canvas", 0, 1, touchstart_callback);
	// emscripten_set_touchend_callback("canvas", 0, 1, touchend_callback);
	// emscripten_set_touchcancel_callback("canvas", 0, 1, touchcancel_callback);
	// emscripten_set_touchmove_callback("canvas", 0, 1, touchmove_callback);


	emscripten_set_main_loop(sdl_init,0,0);
}

}
