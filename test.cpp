#include <cstdio>
#include <cstdlib>
#include "BoltzQVals.hpp"
#include <iostream>

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

int main()
{

    // printf("%20.19e\n",qMaxwell[0][0][0]);
    // printf("%20.19e\n",qHardSphere[0][0][0]);
    // printf("%20.19e\n", computeBoltzQ(10,10,10,false));

    double *y = (double *)std::malloc(N*sizeof(double));
    double *k1 = (double *)std::malloc(N*sizeof(double));
    double *k2 = (double *)std::malloc(N*sizeof(double));
    double *k3 = (double *)std::malloc(N*sizeof(double));
    double *k4 = (double *)std::malloc(N*sizeof(double));
    double *tmp = (double *)std::malloc(N*sizeof(double));

    init(y);
    init(k1);
    init(k2);
    init(k3);
    init(k4);
    init(tmp);

    for(int i = 0; i < N; i++)
    {
        y[i] = 1.0;
    }



    double h = 0.005;
    for(double t = 0; t < 100; t += h)
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

        std::printf("%f\n", t);

        for(int i = 0; i < N; i++)
        {
            std::printf("%20.19e\n",y[i]);
        }

        std::cout << "-----------------" << std::endl;
    }


    std::free(y);
    std::free(k1);
    std::free(k2);
    std::free(k3);
    std::free(k4);
    std::free(tmp);

    return 0;
}
