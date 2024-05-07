#include "dndt.h"
#include "global.h"
#include"E.h"

#include<complex>

std::complex<double> dndt(double t, std::complex<double> p)
{
   std::complex<double> result;

    result = (i / h_cut) * mu * E(t) * (std::conj(p) - p);
    
    return result;
}