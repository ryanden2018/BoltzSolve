#include "BoltzQ.hpp"
#include <cstdio>

int main()
{
    int N = 3;
    printf("#ifndef BOLTZQVALS_HPP\n");
    printf("#define BOLTZQVALS_HPP\n");
    printf("const double boltzQMaxwell[%d][%d][%d] = \n", N, N, N);
    printf("{\n");
    for(int i = 0; i < N; i++)
    {
        printf("    {\n");
        for(int j = 0; j < N; j++)
        {
            printf("        {\n");
            for(int k = 0; k < N; k++)
            {
                printf(k == N-1 ? "            %20.19e\n" : "            %20.19e,\n",computeBoltzQ(i,j,k,false));
            }
            printf(j == N-1 ? "        }\n" : "        },\n");
        }
        printf(i == N-1 ? "    }\n" : "    },\n");
    }
    printf("};\n");
    printf("const double boltzQHardSphere[%d][%d][%d] = \n", N, N, N);
    printf("{\n");
    for(int i = 0; i < N; i++)
    {
        printf("    {\n");
        for(int j = 0; j < N; j++)
        {
            printf("        {\n");
            for(int k = 0; k < N; k++)
            {
                printf(k == N-1 ? "            %20.19e\n" : "            %20.19e,\n",computeBoltzQ(i,j,k,true));
            }
            printf(j == N-1 ? "        }\n" : "        },\n");
        }
        printf(i == N-1 ? "    }\n" : "    },\n");
    }
    printf("};\n");
    printf("#endif\n");
    return 0;
}