#include "GridTrig.hh"
#include "consts.hh"
#include <cmath>

GridTrig::GridTrig(const double s_gp[SDIV+1], const double mu[MDIV+1]) {
    compute_f2n(s_gp);
    compute_f_rho_gamma(s_gp);
    compute_trig(mu);
}

void GridTrig::compute_f2n(const double s_gp[SDIV+1]) {
    for(int n=0; n<=LMAX; ++n) {
        for(int i=2; i<=SDIV; ++i) {
            _f2n[n+1][i] = pow((1.0-s_gp[i])/s_gp[i], 2.0*n);
        }
    }
}

void GridTrig::compute_f_rho_gamma(const double s_gp[SDIV+1]) {
    int j, k, n;
    double sk, sj, sk1, sj1;
    if(SMAX != 1.0) {

        for(j=2; j<=SDIV; ++j) {
            for(n=1; n<=LMAX; ++n) {
                for(k=2; k<=SDIV; ++k) {
                    sk = s_gp[k];
                    sj = s_gp[j];
                    sk1 = 1.0-sk;
                    sj1 = 1.0-sj;

                    if(k < j) {   
                        f_rho[j][n+1][k] = _f2n[n+1][j]*sj1/(sj*_f2n[n+1][k]*sk1*sk1);
                        f_gama[j][n+1][k] = _f2n[n+1][j]/(_f2n[n+1][k]*sk*sk1);
                    } else {     
                        f_rho[j][n+1][k] = _f2n[n+1][k]/(_f2n[n+1][j]*sk*sk1);
                        f_gama[j][n+1][k] = _f2n[n+1][k]*sj1*sj1*sk/(sj*sj*_f2n[n+1][j]*sk1*sk1*sk1);
                    }
                }
            }
        }
        j = 1;
 
        n = 0; 
        for(k=2; k<=SDIV; ++k) {
            sk = s_gp[k];
            f_rho[j][n+1][k] = 1.0/(sk*(1.0-sk));
        }

        n = 1;
        for(k=2; k<=SDIV; ++k) {
            sk = s_gp[k];
            sk1 = 1.0-sk;         
            f_rho[j][n+1][k] = 0.0;
            f_gama[j][n+1][k] = 1.0/(sk*sk1);
        }

        for(n=2; n<=LMAX; ++n) {
            for(k=1; k<=SDIV; ++k) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
            }
        }

        k = 1;
        //Chache inefficient
        n = 0;
        for(j=1; j<=SDIV; ++j) {
            f_rho[j][n+1][k] = 0.0;
        }
        //cache inefficent
        for(j=1; j<=SDIV; ++j) {
            for(n=1; n<=LMAX; ++n) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
            }
        }

        //cache inefficent
        n = 0;
        for(j=2; j<=SDIV; ++j) {
            for(k=2; k<=SDIV; ++k) {
                sk = s_gp[k];
                sj = s_gp[j];
                sk1 = 1.0-sk;
                sj1 = 1.0-sj;

                if(k < j) {
                    f_rho[j][n+1][k] = sj1/(sj*sk1*sk1);
                } else {    
                    f_rho[j][n+1][k] = 1.0/(sk*sk1);
                }
            }
        }         

    } else {      
        for(j=2; j<=SDIV-1; ++j) {
            for(n=1; n<=LMAX; ++n) {
                for(k=2; k<=SDIV-1; ++k) {
                    sk = s_gp[k];
                    sj = s_gp[j];
                    sk1 = 1.0-sk;
                    sj1 = 1.0-sj;

                    if(k < j) {   
                        f_rho[j][n+1][k] = _f2n[n+1][j]*sj1/(sj*_f2n[n+1][k]*sk1*sk1);
                        f_gama[j][n+1][k] = _f2n[n+1][j]/(_f2n[n+1][k]*sk*sk1);
                    } else {     
                        f_rho[j][n+1][k] = _f2n[n+1][k]/(_f2n[n+1][j]*sk*sk1);

                        f_gama[j][n+1][k] = _f2n[n+1][k]*sj1*sj1*sk/(sj*sj*_f2n[n+1][j]*sk1*sk1*sk1);
                    }
                }
            }
        }
   
        j = 1;
 
        n = 0; 
        for(k=2; k<=SDIV-1; ++k) {
            sk = s_gp[k];
            f_rho[j][n+1][k] = 1.0/(sk*(1.0-sk));
        }

        n = 1;
        for(k=2; k<=SDIV-1; ++k) {
            sk = s_gp[k];
            sk1 = 1.0-sk;         
            f_rho[j][n+1][k] = 0.0;
            f_gama[j][n+1][k] = 1.0/(sk*sk1);
        }

        for(n=2; n<=LMAX; ++n) {
            for(k=1; k<=SDIV-1; ++k) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
             }
        }

        k = 1;
 
        //cache inefficient
        n = 0;
        for(j=1; j<=SDIV-1; ++j) {
            f_rho[j][n+1][k]=0.0;
        }

        //cache inefficent
        for(j=1; j<=SDIV-1; ++j) {
            for(n=1; n<=LMAX; ++n) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
            }
        }
 
        //cache inefficient
        n = 0;
        for(j=2; j<=SDIV-1; ++j) {
            for(k=2; k<=SDIV-1; ++k) {
                sk = s_gp[k];
                sj = s_gp[j];
                sk1 = 1.0-sk;
                sj1 = 1.0-sj;

                if(k < j) { 
                  f_rho[j][n+1][k] = sj1/(sj*sk1*sk1);
                } else {     
                  f_rho[j][n+1][k] = 1.0/(sk*sk1);
                }
            }
        }
 
        j = SDIV;
        for(n=1; n<=LMAX; ++n) {
            for(k=1; k<=SDIV; ++k) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
            }
        }

        //cache inefficient
        k = SDIV;
        for(j=1; j<=SDIV; ++j) {
            for(n=1; n<=LMAX; ++n) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
            }
        }
    }
}

void GridTrig::compute_trig(const double mu[MDIV+1]) {
    int n = 0, i, m;
    for(i=1; i<=MDIV; ++i) {
        P_2n[i][n+1] = legendre(2*n, mu[i]);
    }
    
    for(m=1; m<=MDIV; ++m) { 
        sin_theta[m] = sqrt(1.0-mu[m]*mu[m]);  
        _theta[m] = asin(sin_theta[m]);
    }

    for(i=1; i<=MDIV; ++i) {
        for(n=1;n<=LMAX;++n) {
            P_2n[i][n+1] = legendre(2*n, mu[i]);
            P_2n_t[n+1][i] = P_2n[i][n+1];
            P1_2n_1[i][n+1] = plgndr(2*n-1, 1, mu[i]);
            P1_2n_1_t[n+1][i] = P1_2n_1[i][n+1];
            sin_2n_1_theta[i][n] = sin((2.0*n-1.0)*_theta[i]);
            sin_2n_1_theta_t[n][i] = sin_2n_1_theta[i][n];
        }
    }
}
