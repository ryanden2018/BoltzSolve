#include "BoltzQ.hpp"
#include <cstdio>
#include <iostream>

int main()
{
    int N = 11;
    printf("#ifndef BOLTZQVALS_HPP");
    std::cout << std::endl;
    printf("#define BOLTZQVALS_HPP");
    std::cout << std::endl;
    printf("const float boltzQMaxwell[%d][%d][%d] = ", N, N, N);
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
                printf(k == N-1 ? "            %20.19e" : "            %20.19e,",computeBoltzQ(i,j,k,false));
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
    return 0;
}