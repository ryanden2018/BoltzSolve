#include <math.h>
#include <memory.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <SDL/SDL.h>
#include <stdio.h>
#include <stdlib.h>

const int N0 = 20;
const int N = N0*N0;
const double L0 = 16.0;
const double h = L0/N0;
const double PI = 3.141592653589793238462;
const int M = 10;

int clamp(int i)
{
    return i < 0 ? 0 :
        ( i >= N0 ? N0 : i);
}

void computeOp(double *in, double* out)
{
    for(int i = 0; i < N; i++)
    {
        out[i] = 0.0;
        int ix = i%N0;
        int iy = (i-ix)/N0;
        double vx = ix*h-L0/2.0;
        double vy = iy*h-L0/2.0;
        for(int j = 0; j < N; j++)
        {
            int jx = j%N0;
            int jy = (j-jx)/N0;
            double vsx = jx*h-L0/2.0;
            double vsy = jy*h-L0/2.0;
            for(int k = 0; k < M; k++)
            {
                double theta = 2.0 * PI * k / M;
                double tau = cos(theta)*(vx-vsx) + sin(theta)*(vy-vsy);
                double abstau = tau > 0 ? tau : -tau;
                double vpx = vx - cos(theta)*tau;
                double vpy = vy - sin(theta)*tau;
                double vspx = vsx + cos(theta)*tau;
                double vspy = vsy + sin(theta)*tau;
                int ii = clamp(round((vpy+L0/2.0)/h))*N0 + clamp(round((vpx+L0/2.0)/h));
                int jj = clamp(round((vspy+L0/2.0)/h))*N0 + clamp(round((vspx+L0/2.0)/h));
                out[i] += abstau * in[ii] * in[jj] * h*h*2.0*PI/M;
                out[i] -= abstau * in[i] * in[j] * h*h*2.0*PI/M;
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

double sum(double *in)
{
    double val = 0.0;
    for(int i = 0; i < N; i++)
    {
        val += in[i];
    }
    return val;
}

void init(double *out)
{
    for(int i = 0; i < N; i++)
    {
        out[i] = 0.0;
    }
}

void cut(double *out)
{
    for(int i = 0; i < N; i++)
    {
        if(out[i] < 0.0) out[i] = 0.0;
    }
}

double evalData(double x, double y)
{
    double r = sqrt(x*x+y*y);
    //return exp(-r*r) * (1.0-0.75*cos(5*r));
    return (r < 2.0 ? 0.25 : 0.0) - (r<1.0 ? 0.5*0.25 : 0.0);
   // return exp(-r*r) * (1.0+0.75*exp(-4*r*r)*cos(10*r));
    //return exp(-r*r) +0.25* exp(-(r+3)*(r+3)) +0.25* exp(-(r-3)*(r-3));
}

SDL_Surface *screen;

double *y;
double *k1;
double *k2;
double *k3;
double *k4;
double *tmp;
const int NUMPIXELS = 693;

void repaint()
{
	if (SDL_MUSTLOCK(screen)) SDL_LockSurface(screen);
	for (int i = 0; i < NUMPIXELS; i++) {
		for (int j = 0; j < NUMPIXELS; j++) {
			*((Uint32*)screen->pixels + i * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 255, 255, 255, 255);
		}
	}
    for (int j = 0; j < NUMPIXELS; j++)
    {
        double vy = 7.5*(j-0.5*NUMPIXELS)/NUMPIXELS;
        double vx = 0.0;
        int ix = clamp(round((vx+L0/2.0)/h));
        int iy = clamp(round((vy+L0/2.0)/h));
        double val = y[iy*N0+ix];
        int i = (int) (((val+1.0)/5.0)*NUMPIXELS);
        int ii = (int) (((1.0)/5.0)*NUMPIXELS);
        i = NUMPIXELS - i;
        ii = NUMPIXELS - ii;
        if(i < 0) i = 0;
        if(i > NUMPIXELS) i = NUMPIXELS;
        *((Uint32*)screen->pixels + i * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 0, 0, 0, 255);
        *((Uint32*)screen->pixels + ii * NUMPIXELS + j) = SDL_MapRGBA(screen->format, 0, 0, 0, 255);
    }
    printf("%7.6e\n",sum(y));
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

    double dt = 0.05;
    computeOp(y,k1);
    init(tmp);
    add(k1,tmp);
    scale(tmp,dt/2.0);
    add(y,tmp);
    computeOp(tmp,k2);
    init(tmp);
    add(k2,tmp);
    scale(tmp,dt/2.0);
    add(y,tmp);
    computeOp(tmp,k3);
    init(tmp);
    add(k3,tmp);
    scale(tmp,dt);
    add(y,tmp);
    computeOp(tmp,k4);
    init(tmp);
    add(k1,tmp);
    add(k2,tmp);
    add(k2,tmp);
    add(k3,tmp);
    add(k3,tmp);
    add(k4,tmp);
    scale(tmp,dt/6.0);
    add(y,tmp);
    swap(y,tmp);
    cut(y);
}

int main(int argc, char ** argv)
{
    y = (double *)malloc(N*sizeof(double));
    k1 = (double *)malloc(N*sizeof(double));
    k2 = (double *)malloc(N*sizeof(double));
    k3 = (double *)malloc(N*sizeof(double));
    k4 = (double *)malloc(N*sizeof(double));
    tmp = (double *)malloc(N*sizeof(double));

    init(y);
    init(k1);
    init(k2);
    init(k3);
    init(k4);
    init(tmp);

    for(int i = 0; i < N0; i++)
    {
        for(int j = 0; j < N0; j++)
        {
            double vx = j*h-L0/2.0;
            double vy = i*h-L0/2.0;
            y[i*N0+j] = evalData(vx,vy);
        }
    }


    // free(y);
    // free(k1);
    // free(k2);
    // free(k3);
    // free(k4);
    // free(tmp);
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


	emscripten_set_main_loop(sdl_init,10,0);
}
