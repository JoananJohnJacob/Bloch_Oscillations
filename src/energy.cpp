#include"energy.h"
#include<vector>
#include"global.h"

std::vector<double> energy(std::vector<double> &x)
{
    int n = x.size();

     std::vector<double> energy(n);

    
    for(int j=0; j<n; j++)
    {
        energy[j] = h_ev*x[j];
    
    }

    return energy;
}

