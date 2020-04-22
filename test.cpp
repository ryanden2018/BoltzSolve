#include <cstdio>
#include <cstdlib>
#include "BoltzQ.hpp"
int main()
{

    // printf("%20.19e\n",qMaxwell[0][0][0]);
    // printf("%20.19e\n",qHardSphere[0][0][0]);
    printf("%20.19e\n", computeBoltzQ(10,10,10,false));

    return 0;
}
