#include"real_part.h"
#include<cmath>
#include<vector>
#include<iostream>

std::vector<double> real_part(std::vector<std::complex<double>> &x)
{
    
    int n = x.size();
    std::vector<double> real(n);
    
    for(int k = 0; k<n; k++)
    {
        real[k] = x[k].real();
    }

    return real;
}

std::vector<double> imag_part(std::vector<std::complex<double>> &x)
{
    
    int n = x.size();
    std::vector<double> imag(n);

    for(int k = 0; k<n; k++)
    {
        imag[k] = x[k].imag();
    }

    return imag;
}

std::vector<double> abso(std::vector<std::complex<double>> &x)
{
    
    int n = x.size();
    std::vector<double> abso(n);

    for(int k = 0; k<n; k++)
    {
        abso[k] = std::abs(x[k]);
    }

    return abso;
}

std::vector<std::complex<double>> divide(std::vector<std::complex<double>> &x,std::vector<std::complex<double>> &y)
{
    std::vector<std::complex<double>> result(x.size(),0);
    int n1 = x.size(), n2 = y.size();

    if(n1!=n2)
    {
        std::cerr<<"Error: vectors should be of same size"<<std::endl;
        
    
    }
    else
    {
        for(int k=0; k<n1; k++)
        {
            result[k] = x[k]/y[k];
        }
    }

    return result;
}