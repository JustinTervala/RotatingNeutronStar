#pragma once

#include "equil_util.h"

struct GridTrig {
    public:
        double sin_theta[MDIV+1];
        native_tensor<float, SDIV+1, LMAX+2, SDIV+1> f_rho;
        native_tensor<float, SDIV+1, LMAX+2, SDIV+1> f_gama;
        native_matrix<double, MDIV+1, LMAX+2> P_2n; 
        native_matrix<double, LMAX+2, MDIV+1> P_2n_t;
        native_matrix<double, MDIV+1, LMAX+2> P1_2n_1;
        native_matrix<double, LMAX+2, MDIV+1> P1_2n_1_t;
        native_matrix<double, MDIV+1, LMAX+1> sin_2n_1_theta;
        native_matrix<double, LMAX+1, MDIV+1> sin_2n_1_theta_t;

        GridTrig(const double s_gp[SDIV+1], const double mu[MDIV+1]);
        ~GridTrig() {}     
    private:
        native_matrix<double, LMAX+2, SDIV+1>  _f2n;
        double _theta[MDIV+1];
        void compute_f2n(const double s_gp[SDIV+1]);

        void compute_f_rho_gamma(const double s_gp[SDIV+1]);

        void compute_trig(const double mu[MDIV+1]);

};
