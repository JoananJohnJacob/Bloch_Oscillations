#ifndef FOURIER_H
#define FOURIER_H

#include <complex>
#include <vector>

std::vector<std::complex<double>> fft(std::vector<std::complex<double>> &x);
std::vector<std::complex<double>> fft(std::vector<double> &x);
std::vector<std::complex<double>> ifft(std::vector<std::complex<double>> &y);
std::vector<std::complex<double>> ifft(std::vector<double> &y);
std::vector<double> fftfreq(int N, double d=1.0);

#endif
