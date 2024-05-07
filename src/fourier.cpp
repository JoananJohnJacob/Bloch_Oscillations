#include "fourier.h"
#include <fftw3.h>

// Function to compute FFT using FFTW
std::vector<std::complex<double>> fft(std::vector<std::complex<double>> &x) {
    // ... same as before ...
    int N = x.size();
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for(int i = 0; i < N; i++) {
        in[i][0] = x[i].real();
        in[i][1] = x[i].imag();
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    std::vector<std::complex<double>> y(N);
    for(int i = 0; i < N; i++) {
        y[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return y;

}

std::vector<std::complex<double>> fft(std::vector<double> &x) {
    std::vector<std::complex<double>> x_complex(x.begin(), x.end());
    return fft(x_complex);
}

// Function to compute IFFT using FFTW
std::vector<std::complex<double>> ifft(std::vector<std::complex<double>> &y) {
    int N = y.size();
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for(int i = 0; i < N; i++) {
        in[i][0] = y[i].real();
        in[i][1] = y[i].imag();
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    std::vector<std::complex<double>> x(N);
    for(int i = 0; i < N; i++) {
        x[i] = std::complex<double>(out[i][0], out[i][1]) / static_cast<double>(N);
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return x;
}

std::vector<std::complex<double>> ifft(std::vector<double> &y) {
    std::vector<std::complex<double>> y_complex(y.begin(), y.end());
    return ifft(y_complex);
}

// Function to compute FFT frequencies
std::vector<double> fftfreq(int N, double d) {
    double val = 1.0 / (N * d);
    int N_half = N / 2;
    std::vector<double> result(N);
    for(int i = 0; i < N_half; i++) result[i] = i * val;
    for(int i = -N_half; i < 0; i++) result[N+i] = i * val;
    return result;
}
