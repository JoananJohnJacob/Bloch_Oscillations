#include "E.h"
#include "global.h"
#include <cmath>

double E(double t) {
    double E_f, Gaussian, Cosine;
    double t0 = 25e-15, Dt = 5e-15;
    double energy = 2;
    double wL = energy/h_cut;
    double E0 = 100;
    double num =  t-t0;
    Gaussian = E0*exp(-(num*num/(Dt*Dt)));
    Cosine = cos(wL*t);
    E_f = Gaussian*Cosine;
    return E_f;
}

double F(double t){
    double F = (0.2)/(a);

    return F;
}