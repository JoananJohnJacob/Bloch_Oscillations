#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include "matplotlibcpp.h"
#include "dpdt.h"
#include "global.h"
#include "dndt.h"
#include "E.h"
#include "fourier.h"
#include "csv_writer.h"
#include "real_part.h"
#include "energy.h"

namespace plt = matplotlibcpp;

int main() {
    //----------------------------------------------------- INITIALIZING PARAMETERS IN K SPACE --------------------------------------------------------

    //Number of k-points
    const int n_k = 750;
    

    // band Gap
    const double E_gap = 1.5;

    //bandwidth of valence band and conduction band
    const double Delta_E_v = 0.15  , Delta_E_c = 0.85;

    // Declaring vectors for Energies,frequencies and k

    std::vector<double> E_v(n_k); // valence band
    std::vector<double> E_c(n_k); // conduction band 
    std::vector<double> w_k(n_k); // Frequencies
    std::vector<double> kk(n_k); // k vector

    // Initializing parameters 
    
    for (int j = 0; j < n_k; j++){

        //initializing k vector
        kk[j] = (j*(2*pi/n_k) - pi);

        // Energy dispersion
        E_v[j] = Delta_E_v/2*(std::cos(kk[j]) - 1);
        E_c[j] = Delta_E_c/2*(1 -  std::cos(kk[j])) + E_gap;

        //frequency spectrum
        w_k[j] = (E_c[j] - E_v[j]) / h_bar;
    }

    //-------------------------------------------------------------------------------------------------------------------------------------------

    //---------------------------------------------------INITIALIZING TIME AND ELECTRIC FIELD----------------------------------------------------
    // starting and stopping time
    double t_start =0, t_stop = 5000e-15; 
    
    //time steps dt and dt/2 for RK scheme;
    double dt = 0.1e-15, dt2 = dt/2;
    const int N = static_cast<int> ((t_stop - t_start)/dt);

    //Initializing time parameters for the vectors

    double n_t = N;



    
    //declaring time vector and electric field vector

    std::vector<double> t(n_t), E_f(n_t);

    // initializing time and electric field

    for(int j = 0; j < n_t; j++){

        //time
        t[j] = j*dt;
        //electric field
        E_f[j] = E(t[j]);
    }


    /* plt::plot(t,E_f);
    plt::show(); */

    //-------------------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------- INITIALIZING VECTORS ------------------------------------------------------------

    std::vector<std::vector<std::complex<double>>> p; //microscopic polarization vector

    //initializing polarization vector
    p.resize(n_t,std::vector<std::complex<double>>(n_k)); 

    std::vector<std::complex<double>> P_total(n_t); //macroscopic polarization

    //--------------------------------------------------------------------------------------------------------------------------------
    
    //---------------------------------------------------RK-SCHEME--------------------------------------------------------------------
    
    std::complex<double> k1,k2,k3,k4;//values for ks
    std::complex<double> p_now, p_l, p_r, dpdk; // values to store ps 
    double t_now,dk = kk[1]-kk[0];

    for(int time_step = 0; time_step < (n_t-1) ; time_step++){

        for(int k_step = 0; k_step < n_k ; k_step++){

            t_now = t[time_step];
            p_now = p[time_step][k_step];
            //finding p_l and p_r
            if(k_step == 0)
            {
                p_l = p[time_step][n_k];
                p_r = p[time_step][k_step+1];
            }
            else if(k_step == n_k)
            {
                p_l = p[time_step][k_step-1];
                p_r = p[time_step][0];
            }
            else
            {
                p_l =  p[time_step][k_step-1];
                p_r = p[time_step][k_step+1];

            }

            dpdk = (p_r - p_l)/(2*dk);

            k1 = dpdt(t_now, p_now , w_k[k_step],dpdk);
            k2 = dpdt(t_now+dt2, p_now + dt2*k1, w_k[k_step],dpdk);
            k3 = dpdt(t_now+dt2, p_now + dt2*k2, w_k[k_step],dpdk);
            k4 = dpdt(t_now+dt, p_now + dt*k3, w_k[k_step],dpdk);

            p[time_step+1][k_step] = p[time_step][k_step] + (dt/6)*(k1 + 2.0 *k2 + 2.0 *k3 + k4);

        }
    }

    //finding macroscopic polarization
    for(int time_step = 0; time_step < n_t; time_step++){
        for(int k_index = 0 ; k_index < n_k; k_index++){

            P_total[time_step] += 2.0*mu*p[time_step][k_index].real();
        }
    }

    //----------------------------------------------------------FOURIER TRANSFORM --------------------------------------------------------- 

    // frequency spectrum
    std::vector<double> frequency = fftfreq(n_t,dt);
    std::vector<double> Energy = energy(frequency);


    //fourier of total polarization
    std::vector<std::complex<double>> P_w = ifft(P_total);

    //fourier of Electric field
    std::vector<std::complex<double>> E_w = ifft(E_f);

    // Susceptibility

    std::vector<std::complex<double>> alpha_w = divide(P_w,E_w);

    //Absorption
    std::vector<double> alpha_w_i = imag_part(alpha_w);

    plt::figure(1);
    plt::named_plot("Absorption",Energy,alpha_w_i);
    plt::xlabel("Energy in eV");
    plt::ylabel("Absorption");
    plt::legend();
    plt::ylim(0.0,10e-53);
    plt::xlim(1.4,2.6);
    plt::grid("True");

   
    plt::figure(2);
    plt::named_plot("Valence Band",kk,E_v);
    plt::named_plot("Conduction Band",kk,E_c);
    plt::xlabel("k");
    plt::ylabel("E");
    plt::legend();
    //plt::plot(kk,w_k);
    plt::grid("True");
    plt::show();



    return 0;
}
