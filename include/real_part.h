#ifndef REAL_PARTH
#define REAL_PARTH
#include<cmath>
#include<vector>
#include<complex>
#include<iostream>

std::vector<double> real_part(std::vector<std::complex<double>> &x);
std::vector<double> imag_part(std::vector<std::complex<double>> &x);
std::vector<double> abso(std::vector<std::complex<double>> &x);
std::vector<std::complex<double>> divide(std::vector<std::complex<double>> &x,std::vector<std::complex<double>> &y);

#endif

