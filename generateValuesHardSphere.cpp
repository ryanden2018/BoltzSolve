#include "BoltzQ.hpp"
#include <cstdio>
#include <iostream>

int main()
{
    int N = 11;
    printf("const double boltzQHardSphere[%d][%d][%d] = ", N, N, N);
    std::cout << std::endl;
    printf("{");
    std::cout << std::endl;
    for(int i = 0; i < N; i++)
    {
        printf("    {");
        std::cout << std::endl;
        for(int j = 0; j < N; j++)
        {
            printf("        {");
            std::cout << std::endl;
            for(int k = 0; k < N; k++)
            {
                printf(k == N-1 ? "            %20.19e" : "            %20.19e,",computeBoltzQ(i,j,k,true));
                std::cout << std::endl;
            }
            printf(j == N-1 ? "        }" : "        },");
            std::cout << std::endl;
        }
        printf(i == N-1 ? "    }" : "    },");
        std::cout << std::endl;
    }
    printf("};");
    std::cout << std::endl;
    printf("#endif");
    std::cout << std::endl;
    return 0;
}