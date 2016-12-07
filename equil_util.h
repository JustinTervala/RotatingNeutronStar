#include <array>
#include "consts.h"
template <class T, size_t size_x, size_t size_y>
std::array<std::array<T, size_y>, size_x> copy_from_nr(T** in) {
    std::array<std::array<T, size_y>, size_x> tmp;
    for (size_t ii = 0; ii < size_x; ++ii) {
        for (size_t jj = 0; jj < size_y; ++jj) {
            tmp[ii][jj] = in[ii+1][jj+1];    
        }
    }
    return tmp;
}

template<class T, size_t size_x, size_t size_y>
void copy_to_nr(std::array<std::array<T, size_y>, size_x> in, T** out) {
    for (size_t ii = 0; ii < size_x; ++ii) {
        for (size_t jj = 0; jj < size_y; ++jj) {
            out[ii+1][jj+1] = in[ii][jj];    
        }
    }
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

double legendre(int n, double x);

double plgndr(int l, int m, double x);
 
double rtsec_G(double (*func)(double, double), 
               double Gamma_P,
               double x1, 
               double x2, 
               double xacc,
               double ee);

