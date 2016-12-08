#include <array>
#include <string>
#include "consts.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

struct Metric {
    double rho, gama, alpha, omega;
    Metric() : rho(0.0), gama(0.0), alpha(0.0), omega(0.0) {}
    ~Metric() {}
};

template <typename T> int sgn(T val, T& out) {
    out = (T(0) < val) - (val < T(0));
}

void hunt(double xx[], int n, double x, int *jlo);

double interp(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt);

double deriv_s(double **f, int s, int m);

template <class T, size_t size_x, size_t size_y>
T deriv_s(const std::array<std::array<T, size_y>, size_x>& f, size_t s, size_t m) {
    double d_temp;

    switch(s) { 
        case 1: d_temp = (f[s+1][m]-f[s][m])/DS;
                break;

        case SDIV: d_temp = (f[s][m]-f[s-1][m])/DS;
                   break;
      
        default: d_temp = (f[s+1][m]-f[s-1][m])/(2.0*DS);
                 break; 
    } 
    return d_temp;
}

template <size_t size_x, size_t size_y>
double deriv_s(const std::array<std::array<Metric, size_y>, size_x>& metric, double Metric::*field, size_t s, size_t m){
    double d_temp;

    switch(s) { 
        case 1: d_temp = (metric[s+1][m].*field-metric[s][m].*field)/DS;
                break;

        case SDIV: d_temp = (metric[s][m].*field-metric[s-1][m].*field)/DS;
                   break;
      
        default: d_temp = (metric[s+1][m].*field-metric[s-1][m].*field)/(2.0*DS);
                 break; 
    } 
    return d_temp;
}

double deriv_ss(double **f, int s, int m);

double deriv_m(double **f, int s, int m);

template <class T, size_t size_x, size_t size_y>
T deriv_m(const std::array<std::array<T, size_y>, size_x>& f, size_t s, size_t m) {
    double d_temp;

    switch(m) { 
        case 1: d_temp = (f[s][m+1]-f[s][m])/DM;
                break; 

        case MDIV: d_temp = (f[s][m]-f[s][m-1])/DM;
                   break;
      
        default: d_temp = (f[s][m+1]-f[s][m-1])/(2.0*DM);
                 break; 
    } 
    return d_temp;
}

template <size_t size_x, size_t size_y>
double deriv_m(const std::array<std::array<Metric, size_y>, size_x>& metric, double Metric::*field, size_t s, size_t m) {
    double d_temp;

    switch(m) { 
        case 1: d_temp = (metric[s][m+1].*field-metric[s][m].*field)/DM;
                break; 

        case MDIV: d_temp = (metric[s][m].*field-metric[s][m-1].*field)/DM;
                   break;
      
        default: d_temp = (metric[s][m+1].*field-metric[s][m-1].*field)/(2.0*DM);
                 break; 
    } 
    return d_temp;
}

double deriv_mm(double **f, int s, int m);

double deriv_sm(double **f, int s, int m);

template<class T, size_t size_x, size_t size_y>
T deriv_sm(const std::array<std::array<T, size_y>, size_x>& f, int s, int m) {
    double d_temp;

    switch(s) {
        case 1: 
            if(m == 1) {   
                d_temp = (f[s+1][m+1]-f[s][m+1]-f[s+1][m]+f[s][m])/(DM*DS);
            } else {
                if(m == MDIV) {
                    d_temp = (f[s+1][m]-f[s][m]-f[s+1][m-1]+f[s][m-1])/(DM*DS);
                } else {         
                    d_temp = (f[s+1][m+1]-f[s+1][m-1]-f[s][m+1]+f[s][m-1])/(2.0*DM*DS);
                }
            }
            break;

        case SDIV: 
            if(m == 1) {   
                d_temp = (f[s][m+1]-f[s][m]-f[s-1][m+1]+f[s-1][m])/(DM*DS);
            } else {
                if(m == MDIV) {
                    d_temp = (f[s][m]-f[s-1][m]-f[s][m-1]+f[s-1][m-1])/(DM*DS);
                } else {         
                    d_temp = (f[s][m+1]-f[s][m-1]-f[s-1][m+1]+f[s-1][m-1])/(2.0*DM*DS);
                }
            }
            break;
  
        default: 
            if(m == 1) {   
                d_temp = (f[s+1][m+1]-f[s-1][m+1]-f[s+1][m]+f[s-1][m])/(2.0*DM*DS);
            } else {
                if(m == MDIV) {
                    d_temp = (f[s+1][m]-f[s-1][m]-f[s+1][m-1]+f[s-1][m-1])/(2.0*DM*DS);
                } else {         
                  d_temp = (f[s+1][m+1]-f[s-1][m+1]-f[s+1][m-1]+f[s-1][m-1])/(4.0*DM*DS);
                }
            }
            break;
    }

    return d_temp;
}

template<size_t size_x, size_t size_y>
double deriv_sm(const std::array<std::array<Metric, size_y>, size_x>& metric, double Metric::*field, int s, int m) {
    double d_temp;

    switch(s) {
        case 1: 
            if(m == 1) {   
                d_temp = (metric[s+1][m+1].*field-metric[s][m+1].*field-metric[s+1][m].*field+metric[s][m].*field)/(DM*DS);
            } else {
                if(m == MDIV) {
                    d_temp = (metric[s+1][m].*field-metric[s][m].*field-metric[s+1][m-1].*field+metric[s][m-1].*field)/(DM*DS);
                } else {         
                    d_temp = (metric[s+1][m+1].*field-metric[s+1][m-1].*field-metric[s][m+1].*field+metric[s][m-1])/(2.0*DM*DS);
                }
            }
            break;

        case SDIV: 
            if(m == 1) {   
                d_temp = (metric[s][m+1].*field-metric[s][m].*field-metric[s-1][m+1].*field+metric[s-1][m].*field)/(DM*DS);
            } else {
                if(m == MDIV) {
                    d_temp = (metric[s][m].*field-metric[s-1][m].*field-metric[s][m-1].*field+metric[s-1][m-1].*field)/(DM*DS);
                } else {         
                    d_temp = (metric[s][m+1].*field-metric[s][m-1].*field-metric[s-1][m+1].*field+metric[s-1][m-1].*field)/(2.0*DM*DS);
                }
            }
            break;
  
        default: 
            if(m == 1) {   
                d_temp = (metric[s+1][m+1].*field-metric[s-1][m+1].*field-metric[s+1][m].*field+metric[s-1][m].*field)/(2.0*DM*DS);
            } else {
                if(m == MDIV) {
                    d_temp = (metric[s+1][m].*field-metric[s-1][m].*field-metric[s+1][m-1].*field+metric[s-1][m-1].*field)/(2.0*DM*DS);
                } else {         
                  d_temp = (metric[s+1][m+1].*field-metric[s-1][m+1].*field-metric[s+1][m-1].*field+metric[s-1][m-1].*field)/(4.0*DM*DS);
                }
            }
            break;
    }

    return d_temp;
}

double legendre(int n, double x);

double plgndr(int l, int m, double x);
 
double rtsec_G(double (*func)(double, double), 
               double Gamma_P,
               double x1, 
               double x2, 
               double xacc,
               double ee);

void print_error(const std::string& error);
