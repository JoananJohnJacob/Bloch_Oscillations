#include "dpdt.h"
#include "global.h"
#include<cmath>
#include "E.h"

std::complex<double> dpdt(double t, std::complex<double> p, double omega, std::complex<double> dpdk){
    std::complex<double> result;
    //double omega = Energy / h_cut;

    result = -i*omega*p + (i/h_bar)*mu*E(t) - (p/T2) + (i*e*F(t)*dpdk)/(i*h_bar);


    return result;

}