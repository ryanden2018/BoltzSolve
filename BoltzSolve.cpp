#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <Eigen/IterativeLinearSolvers>
#include <stdio.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <SDL/SDL.h>
#include <queue>
#include <unsupported/Eigen/FFT>

#include "GaussQuad.hpp"
#include "BoltzQ.hpp"
#include "Hermite.hpp"

const double PI = 3.141592653589793238462;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Triplet<double> Trip;
typedef Eigen::MatrixXd Mat;

// precomputable quantities
const int POLYMAX = 10;
double weights[22];
double coords[22];
double normLegendreDerivProducts[POLYMAX+1][POLYMAX+1];
double normLegendreAltProducts[POLYMAX+1][POLYMAX+1];
double normLegendreLeftVals[POLYMAX+1];
double normLegendreRightVals[POLYMAX+1];
double normLegendreDerivLeftVals[POLYMAX+1];
double normLegendreDerivRightVals[POLYMAX+1];

// double LegendreEval(int p, double y)
// {
// 	if(p == 0) return 1.0;
// 	if(p == 1) return y;
// 	double prev = 1.0;
// 	double cur = y;
// 	for(int n = 1; n < p; n++)
// 	{
// 		double next = ((2.0*n+1.0)*y*cur - n*prev)/(n+1.0);
// 		prev = cur;
// 		cur = next;
// 	}
// 	return cur;
// }

// double LegendreDerivEval(int p, double y)
// {
// 	double val = 0.0;
// 	for(int n = 0; n < p; n++)
// 	{
// 		val = (n+1.0)*LegendreEval(n,y) + y*val;
// 	}
// 	return val;
// }

// double LegendreL2Norm(int p)
// {
// 	return std::pow(2.0/(2.0*p+1.0),0.5);
// }



// double LegendreEvalNorm(int p, double y)
// {
// 	return LegendreEval(p,y) / LegendreL2Norm(p);
// }

// double LegendreDerivEvalNorm(int p, double y)
// {
// 	return LegendreDerivEval(p,y) / LegendreL2Norm(p);
// }

// void MakeLegendreEndpointVals()
// {
// 	for(int p = 0; p < POLYMAX+1; p++)
// 	{
// 		normLegendreLeftVals[p] = LegendreEvalNorm(p,-1.0);
// 		normLegendreRightVals[p] = LegendreEvalNorm(p,1.0);
// 		normLegendreDerivLeftVals[p] = LegendreDerivEvalNorm(p,-1.0);
// 		normLegendreDerivRightVals[p] = LegendreDerivEvalNorm(p,1.0);
// 	}
// }

// void MakeLegendreDerivProducts()
// {
// 	for(int p = 0; p < POLYMAX+1; p++)
// 	{
// 		for(int q = 0; q < POLYMAX+1; q++)
// 		{
// 			double res = 0.0;
// 			for(int k = 0; k < 22; k++)
// 			{
// 				res += weights[k] * LegendreDerivEvalNorm(p,coords[k]) * LegendreDerivEvalNorm(q,coords[k]);
// 			}
// 			normLegendreDerivProducts[p][q] = res;
// 		}
// 	}
// }

// void MakeLegendreAltProducts()
// {
// 	for(int p = 0; p < POLYMAX+1; p++)
// 	{
// 		for(int q = 0; q < POLYMAX+1; q++)
// 		{
// 			double res = 0.0;
// 			for(int k = 0; k < 22; k++)
// 			{
// 				res += weights[k] * LegendreEvalNorm(p,coords[k]) * LegendreDerivEvalNorm(q,coords[k]);
// 			}
// 			normLegendreAltProducts[p][q] = res;
// 		}
// 	}
// }


extern "C" {

const int NUMPIXELS = 693;

double len = 1.0;
SDL_Surface *screen;


void mouse_update(const EmscriptenMouseEvent *e)
{
}

void touch_update(const EmscriptenTouchEvent *e)
{
}

const int NUMLEVELS = 13;

int getColorIndex(double val)
{
	return std::floor( (NUMLEVELS-1) * ((val+0.0001)/1.0002) );
}

double getLambda(double val, int colorIndex)
{
	return (val-colorIndex/(1.0*NUMLEVELS-1.0))*(1.0*NUMLEVELS-1.0);
}

double red(double lambda, int colorIndex, double low = 1.0, double high = 254.0)
{
	return 0.47<lambda && lambda<0.53 ? high : low;
}

double green(double lambda, int colorIndex, double low = 1.0, double high = 254.0)
{
	return 0.47<lambda && lambda<0.53 ? high : low;
}

double blue(double lambda, int colorIndex, double low = 1.0, double high = 254.0)
{
	return 0.47<lambda && lambda<0.53 ? high : low;
}

// void repaintHigh()
// {

// 	if (SDL_MUSTLOCK(screen)) SDL_LockSurface(screen);
// 	for (int i = 0; i < NUMPIXELS; i++) {
// 		for (int j = 0; j < NUMPIXELS; j++) {
// 			double val = (convDiffHigh.Eval(len*(1.0*j)/NUMPIXELS,len*(1.0*i)/NUMPIXELS)-minphi)/(maxphi-minphi);
// 			val = 0.97*(val-0.5)+0.5;
// 			int colorIndex = getColorIndex(val);
// 			double lambda = getLambda(val,colorIndex);
// 			double valr = red(lambda,colorIndex)*255.0;
// 			double valg = green(lambda,colorIndex)*255.0;
// 			double valb = blue(lambda,colorIndex)*255.0;
// 			*((Uint32*)screen->pixels + i * NUMPIXELS + j) = SDL_MapRGBA(screen->format, (int)valr, (int)valg, (int)valb, 255);
// 		}
// 	}
// 	if (SDL_MUSTLOCK(screen)) SDL_UnlockSurface(screen);
// 	SDL_Flip(screen); 
// }


// void repaintLow()
// {

// 	for(int i = 0; i < 11; i++)
// 	{
// 		for(int j = 0; j < 11; j++)
// 		{
// 			dispTemp(i,j) = (convDiff.Eval(len*(1.0*(j+0.5))/11,len*(1.0*(i+0.5))/11));
// 		}
// 	}

// 	FFT2D(dispTemp,dispTempFFTRe,dispTempFFTIm);

// 	for(int i = 0; i < 6; i++)
// 	{
// 		for(int j = 0; j < 6; j++)
// 		{
// 			dispFFTRe(i,j) = dispTempFFTRe(i,j);
// 			dispFFTIm(i,j) = dispTempFFTIm(i,j);
// 		}
// 	}

// 	for(int i = 1; i < 6; i++)
// 	{
// 		for(int j = 0; j < 6; j++)
// 		{
// 			dispFFTRe(NUMPIXELS-i,j) = dispTempFFTRe(11-i,j);
// 			dispFFTIm(NUMPIXELS-i,j) = dispTempFFTIm(11-i,j);
// 		}
// 	}

// 	for(int i = 0; i < 6; i++)
// 	{
// 		for(int j = 1; j < 6; j++)
// 		{
// 			dispFFTRe(i,NUMPIXELS-j) = dispTempFFTRe(i,11-j);
// 			dispFFTIm(i,NUMPIXELS-j) = dispTempFFTIm(i,11-j);
// 		}
// 	}

// 	for(int i = 1; i < 6; i++)
// 	{
// 		for(int j = 1; j < 6; j++)
// 		{
// 			dispFFTRe(NUMPIXELS-i,NUMPIXELS-j) = dispTempFFTRe(11-i,11-j);
// 			dispFFTIm(NUMPIXELS-i,NUMPIXELS-j) = dispTempFFTIm(11-i,11-j);
// 		}
// 	}

// 	IFFT2D(dispFFTRe,dispFFTIm,dispRe,dispIm);

// 	double maxphi = dispRe.maxCoeff();
// 	double minphi = dispRe.minCoeff();

	

// 	if (SDL_MUSTLOCK(screen)) SDL_LockSurface(screen);
// 	for (int i = 0; i < NUMPIXELS; i++) {
// 		for (int j = 0; j < NUMPIXELS; j++) {
// 			double val = (dispRe((i+NUMPIXELS-32)%NUMPIXELS,(j+NUMPIXELS-32)%NUMPIXELS)-minphi)/(maxphi-minphi);
// 			val = 0.97*(val-0.5)+0.5;
// 			int colorIndex = getColorIndex(val);
// 			double lambda = getLambda(val,colorIndex);
// 			double valr = red(lambda,colorIndex,1.0,175.0)*255.0;
// 			double valg = green(lambda,colorIndex,1.0,175.0)*255.0;
// 			double valb = blue(lambda,colorIndex,1.0,175.0)*255.0;
// 			*((Uint32*)screen->pixels + i * NUMPIXELS + j) = SDL_MapRGBA(screen->format, (int)valr, (int)valg, (int)valb, 255);
// 		}
// 	}
// 	if (SDL_MUSTLOCK(screen)) SDL_UnlockSurface(screen);
// 	SDL_Flip(screen); 
// }


EM_BOOL mouseclick_callback(int eventType, const EmscriptenMouseEvent *e, void *userData)
{
	mouse_update(e);
	return 1;
}

EM_BOOL mouseleave_callback(int eventType, const EmscriptenMouseEvent *e, void *userData)
{
	return 1;
}

EM_BOOL mouseup_callback(int eventType, const EmscriptenMouseEvent *e, void *userData)
{
	return 1;
}

EM_BOOL mousedown_callback(int eventType, const EmscriptenMouseEvent *e, void *userData)
{
	mouse_update(e);
	return 1;
}

EM_BOOL mousemove_callback(int eventType, const EmscriptenMouseEvent *e, void *userData)
{
	mouse_update(e);
	return 1;
}

EM_BOOL touchmove_callback(int eventType, const EmscriptenTouchEvent *e, void *userData)
{
	touch_update(e);
	return 1;
}

EM_BOOL touchstart_callback(int eventType, const EmscriptenTouchEvent *e, void *userData)
{
	touch_update(e);
	return 1;
}

EM_BOOL touchend_callback(int eventType, const EmscriptenTouchEvent *e, void *userData)
{
	touch_update(e);
	return 1;
}

EM_BOOL touchcancel_callback(int eventType, const EmscriptenTouchEvent *e, void *userData)
{
	return 1;
}



void init()
{
	
}


int main(int argc, char ** argv)
{
	// MakeWeights();
	// MakeLegendreDerivProducts();
	// MakeLegendreAltProducts();
	// MakeLegendreEndpointVals();

	SDL_Init(SDL_INIT_VIDEO);
	screen = SDL_SetVideoMode(NUMPIXELS, NUMPIXELS, 32, SDL_SWSURFACE);

	emscripten_set_click_callback("canvas", 0, 1, mouseclick_callback);
	emscripten_set_mousedown_callback("canvas", 0, 1, mousedown_callback);
	emscripten_set_mouseup_callback("canvas", 0, 1, mouseup_callback);
	emscripten_set_mousemove_callback("canvas", 0, 1, mousemove_callback);
	emscripten_set_mouseleave_callback("canvas", 0, 1, mouseleave_callback);
	emscripten_set_touchstart_callback("canvas", 0, 1, touchstart_callback);
	emscripten_set_touchend_callback("canvas", 0, 1, touchend_callback);
	emscripten_set_touchcancel_callback("canvas", 0, 1, touchcancel_callback);
	emscripten_set_touchmove_callback("canvas", 0, 1, touchmove_callback);


	emscripten_set_main_loop(init,0,0);
}

}