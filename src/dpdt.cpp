#include "dpdt.h"
#include "global.h"
#include<cmath>
#include "E.h"

std::complex<double> dpdt(double t, std::complex<double> p, double omega, std::complex<double> dpdk){
    std::complex<double> result;
    //double omega = Energy / h_cut;

    result = -i*omega*p + (i/h_cut)*mu*E(t) - gamma_p*p + i*e*F(t)*dpdk;

    return result;

}