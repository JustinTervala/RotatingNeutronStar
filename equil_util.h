#pragma once

#include <array>
#include <string>
#include "consts.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
/*
template<class T, size_t size, size_t... sizes>
struct NdArrayImpl {
    using type = std::array<typename ArrayImpl<T, sizes...>::type, size>;
};

template<class T, size_t size>
struct NdArrayImpl<T, size> {
    using type = std::array<T, size>;
};

template<class T, size_t...sizes>
using ndarray = typename NdArrayImpl<T, sizes...>::type;
*/

template<class T, size_t size_x, size_t size_y>
using matrix = std::array<std::array<T, size_y>, size_x>;
//using matrix = ndarray<T, size_x, size_y>;

template<class T, size_t size_x, size_t size_y>
using native_matrix = T[size_x][size_y];

template<class T, size_t size_x, size_t size_y, size_t size_z>
using tensor = std::array<std::array<std::array<T, size_z>, size_y>, size_x>;

template<class T, size_t size_x, size_t size_y, size_t size_z>
using native_tensor = T[size_x][size_y][size_z];

struct Metric {
    double rho, gama, omega, alpha;
    Metric() : rho(0.0), gama(0.0), omega(0.0), alpha(0.0) {}
    ~Metric() {}
};

struct RhoGamaOmega {
    double rho, gama, omega;
    RhoGamaOmega() : rho(0.0), gama(0.0), omega(0.0) {}
    ~RhoGamaOmega() {}
};

struct ControlConsts {
    const double accuracy;
    const double cf;
    constexpr ControlConsts() : accuracy(0.0), cf(1) {}
    constexpr ControlConsts(double _accuracy, double _cf) : accuracy(_accuracy), cf(_cf) {}
    ~ControlConsts() {} 
};


template <typename T> int sgn(T val, T& out) {
    out = (T(0) < val) - (val < T(0));
}

void hunt(const double xx[], int n, double x, int &jlo);

double interp(const double xp[], 
              const double yp[], 
              int    np ,
              double xb, 
              int    &n_nearest_pt);

double deriv_s(double **f, int s, int m);

template <class T, size_t size_x, size_t size_y>
T deriv_s(const std::array<std::array<T, size_y>, size_x>& f, size_t s, size_t m) {
    switch(s) { 
        case 1: return (f[s+1][m]-f[s][m])/DS;
                break;

        case SDIV: return (f[s][m]-f[s-1][m])/DS;
                   break;
      
        default: return (f[s+1][m]-f[s-1][m])/(2.0*DS);
                 break; 
    } 
}

template <size_t size_x, size_t size_y>
double deriv_s(const std::array<std::array<Metric, size_y>, size_x>& metric, double Metric::*field, size_t s, size_t m){
    switch(s) { 
        case 1: return (metric[s+1][m].*field-metric[s][m].*field)/DS;
                break;

        case SDIV: return (metric[s][m].*field-metric[s-1][m].*field)/DS;
                   break;
      
        default: return (metric[s+1][m].*field-metric[s-1][m].*field)/(2.0*DS);
                 break; 
    } 
}

double deriv_ss(double **f, int s, int m);

double deriv_m(double **f, int s, int m);

template <class T, size_t size_x, size_t size_y>
T deriv_m(const std::array<std::array<T, size_y>, size_x>& f, size_t s, size_t m) {
    switch(m) { 
        case 1: return (f[s][m+1]-f[s][m])/DM;
                break; 

        case MDIV: return (f[s][m]-f[s][m-1])/DM;
                   break;
      
        default: return (f[s][m+1]-f[s][m-1])/(2.0*DM);
                 break; 
    } 
}

template <size_t size_x, size_t size_y>
double deriv_m(const std::array<std::array<Metric, size_y>, size_x>& metric, double Metric::*field, size_t s, size_t m) {
    switch(m) { 
        case 1: return (metric[s][m+1].*field-metric[s][m].*field)/DM;
                break; 

        case MDIV: return (metric[s][m].*field-metric[s][m-1].*field)/DM;
                   break;
      
        default: return (metric[s][m+1].*field-metric[s][m-1].*field)/(2.0*DM);
                 break; 
    } 
}

double deriv_mm(double **f, int s, int m);

double deriv_sm(double **f, int s, int m);

template<class T, size_t size_x, size_t size_y>
T deriv_sm(const std::array<std::array<T, size_y>, size_x>& f, int s, int m) {
    switch(s) {
        case 1: 
            if(m == 1) {   
                return (f[s+1][m+1]-f[s][m+1]-f[s+1][m]+f[s][m])/(DM*DS);
            } else {
                if(m == MDIV) {
                    return (f[s+1][m]-f[s][m]-f[s+1][m-1]+f[s][m-1])/(DM*DS);
                } else {         
                    return (f[s+1][m+1]-f[s+1][m-1]-f[s][m+1]+f[s][m-1])/(2.0*DM*DS);
                }
            }
            break;

        case SDIV: 
            if(m == 1) {   
                return (f[s][m+1]-f[s][m]-f[s-1][m+1]+f[s-1][m])/(DM*DS);
            } else {
                if(m == MDIV) {
                    return (f[s][m]-f[s-1][m]-f[s][m-1]+f[s-1][m-1])/(DM*DS);
                } else {         
                    return (f[s][m+1]-f[s][m-1]-f[s-1][m+1]+f[s-1][m-1])/(2.0*DM*DS);
                }
            }
            break;
  
        default: 
            if(m == 1) {   
                return (f[s+1][m+1]-f[s-1][m+1]-f[s+1][m]+f[s-1][m])/(2.0*DM*DS);
            } else {
                if(m == MDIV) {
                    return (f[s+1][m]-f[s-1][m]-f[s+1][m-1]+f[s-1][m-1])/(2.0*DM*DS);
                } else {         
                  return (f[s+1][m+1]-f[s-1][m+1]-f[s+1][m-1]+f[s-1][m-1])/(4.0*DM*DS);
                }
            }
            break;
    }
}

template<size_t size_x, size_t size_y>
double deriv_sm(const std::array<std::array<Metric, size_y>, size_x>& metric, double Metric::*field, int s, int m) {
    switch(s) {
        case 1: 
            if(m == 1) {   
                return (metric[s+1][m+1].*field-metric[s][m+1].*field-metric[s+1][m].*field+metric[s][m].*field)/(DM*DS);
            } else {
                if(m == MDIV) {
                    return (metric[s+1][m].*field-metric[s][m].*field-metric[s+1][m-1].*field+metric[s][m-1].*field)/(DM*DS);
                } else {         
                    return (metric[s+1][m+1].*field-metric[s+1][m-1].*field-metric[s][m+1].*field+metric[s][m-1].*field)/(2.0*DM*DS);
                }
            }
            break;

        case SDIV: 
            if(m == 1) {   
                return (metric[s][m+1].*field-metric[s][m].*field-metric[s-1][m+1].*field+metric[s-1][m].*field)/(DM*DS);
            } else {
                if(m == MDIV) {
                    return (metric[s][m].*field-metric[s-1][m].*field-metric[s][m-1].*field+metric[s-1][m-1].*field)/(DM*DS);
                } else {         
                    return (metric[s][m+1].*field-metric[s][m-1].*field-metric[s-1][m+1].*field+metric[s-1][m-1].*field)/(2.0*DM*DS);
                }
            }
            break;
  
        default: 
            if(m == 1) {   
                return (metric[s+1][m+1].*field-metric[s-1][m+1].*field-metric[s+1][m].*field+metric[s-1][m].*field)/(2.0*DM*DS);
            } else {
                if(m == MDIV) {
                    return (metric[s+1][m].*field-metric[s-1][m].*field-metric[s+1][m-1].*field+metric[s-1][m-1].*field)/(2.0*DM*DS);
                } else {         
                  return (metric[s+1][m+1].*field-metric[s-1][m+1].*field-metric[s+1][m-1].*field+metric[s-1][m-1].*field)/(4.0*DM*DS);
                }
            }
            break;
    }
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
