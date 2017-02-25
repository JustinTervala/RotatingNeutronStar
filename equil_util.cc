#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <iostream>
#include "consts.hh"
#include "equil_util.hh"
#include <time.h>

long getElapsedTimeNs(struct timespec start, struct timespec stop) {
    return (stop.tv_sec - start.tv_sec)*1e9 + stop.tv_nsec - start.tv_nsec;
}


using namespace consts;


/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/* Adapted from Numerical Recipes.                                         */
/***************************************************************************/
void hunt(const double xx[], int n, double x, int &jlo) { 
    unsigned jm, jhi, inc;
    const bool ascnd = (xx[n] > xx[1]);
    if(jlo <= 0 || jlo > n) {
        jlo = 0;
        jhi = n+1;
    } else {
        if(x >= xx[jlo] == ascnd) {
            if (jlo == n) {
                return;
            }
            inc = 1;
            jhi = jlo+1;
            while(x >= xx[jhi] == ascnd) {
                jlo = jhi;
                inc += inc;
                jhi = (jlo)+inc;
                if (jhi > n) {
                    jhi = n+1;
                    break;
                }
            }
        } else {
            if (jlo == 1) {
                jlo = 0;
                return;
            }
            inc = 1;
            jhi = jlo;
            jlo -= 1;
            while(x < xx[jlo] == ascnd) {
                jhi = jlo;
                inc += inc;
                jlo = jhi-inc;
                if (jlo < 1) {
                    jlo = 0;
                    break;
                }
            }
        }
    }
    while(jhi-jlo != 1) {
        jm = (jhi+jlo) >> 1;
        if (x > xx[jm] == ascnd) {
            jlo = jm;
        } else {
            jhi = jm;
        }
    }
}

/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp(const double xp[], 
              const double yp[], 
              int    np ,
              double xb, 
              int    &n_nearest_pt) { 
    int k,        /* index of 1st point */
        m = 4;      /* degree of interpolation */ 
 
    hunt(xp, np, xb, n_nearest_pt);

    k = std::min(std::max(n_nearest_pt-(m-1)/2, 1), np+1-m);

    if(xb == xp[k] || xb == xp[k+1] || xb == xp[k+2] || xb == xp[k+3]) {
        xb += math::dbl_epsilon;
    }
    double tmp = xb-xp[k];
    
    double y = (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
        +tmp*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
         ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
        +tmp*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
         ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
        +tmp*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
         ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

    return y;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_s(double **f, int s, int m) { 
    double d_temp;

    switch(s) { 
        case 1: d_temp = (f[s+1][m]-f[s][m])/grid::ds;
                break;

        case SDIV: d_temp = (f[s][m]-f[s-1][m])/grid::ds;
                   break;
      
        default: d_temp = (f[s+1][m]-f[s-1][m])/(2.0*grid::ds);
                 break; 
    } 
    return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_ss(double **f, int s, int m) { 
    double d_temp;

    switch(s) { 
        case 1: s = 4;
                d_temp = (f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*square(grid::ds));
                break;

           case 2: s = 4;
                   d_temp = (f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*square(grid::ds));
                   break;

           case 3: s = 4;
                   d_temp = (f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*square(grid::ds));
                   break;

           case SDIV-1: s = SDIV-2;
                        d_temp = (f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*square(grid::ds));
                        break;

           case SDIV: s = SDIV-2;
                      d_temp = (f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*square(grid::ds));
                      break;
      
           default: d_temp = (f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*square(grid::ds));
                    break; 


    }   
    return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. mu                                */ 
/*******************************************************************/
double deriv_m(double **f, int s, int m) {
    double d_temp;

    switch(m) { 
        case 1: d_temp = (f[s][m+1]-f[s][m])/grid::dm;
                break; 

        case MDIV: d_temp = (f[s][m]-f[s][m-1])/grid::dm;
                   break;
      
        default: d_temp = (f[s][m+1]-f[s][m-1])/(2.0*grid::dm);
                 break; 
    } 
    return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_mm(double **f, int s, int m) { 
    double d_temp;

    switch(m) { 
        case 1: m = 2;
                d_temp = (f[s][m+1]-2.0*f[s][m]+f[s][m-1])/square(grid::dm);
                break;

        case MDIV: m = MDIV-1;
                   d_temp = (f[s][m+1]-2.0*f[s][m]+f[s][m-1])/square(grid::dm);
                   break;

        default: d_temp = (f[s][m+1]-2.0*f[s][m]+f[s][m-1])/square(grid::dm);
                 break; 
    } 
    return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s and mu                          */ 
/*******************************************************************/
double deriv_sm(double **f, int s, int m) {
    double d_temp;

    switch(s) {
        case 1: 
            if(m == 1) {   
                d_temp = (f[s+1][m+1]-f[s][m+1]-f[s+1][m]+f[s][m])/(grid::dm*grid::ds);
            } else {
                if(m == MDIV) {
                    d_temp = (f[s+1][m]-f[s][m]-f[s+1][m-1]+f[s][m-1])/(grid::dm*grid::ds);
                } else {         
                    d_temp = (f[s+1][m+1]-f[s+1][m-1]-f[s][m+1]+f[s][m-1])/(2.0*grid::dm*grid::ds);
                }
            }
            break;

        case SDIV: 
            if(m == 1) {   
                d_temp = (f[s][m+1]-f[s][m]-f[s-1][m+1]+f[s-1][m])/(grid::dm*grid::ds);
            } else {
                if(m == MDIV) {
                    d_temp = (f[s][m]-f[s-1][m]-f[s][m-1]+f[s-1][m-1])/(grid::dm*grid::ds);
                } else {         
                    d_temp = (f[s][m+1]-f[s][m-1]-f[s-1][m+1]+f[s-1][m-1])/(2.0*grid::dm*grid::ds);
                }
            }
            break;
  
        default: 
            if(m == 1) {   
                d_temp = (f[s+1][m+1]-f[s-1][m+1]-f[s+1][m]+f[s-1][m])/(2.0*grid::dm*grid::ds);
            } else {
                if(m == MDIV) {
                    d_temp = (f[s+1][m]-f[s-1][m]-f[s+1][m-1]+f[s-1][m-1])/(2.0*grid::dm*grid::ds);
                } else {         
                  d_temp = (f[s+1][m+1]-f[s-1][m+1]-f[s+1][m-1]+f[s-1][m-1])/(4.0*grid::dm*grid::ds);
                }
            }
            break;
    }

    return d_temp;

}

/*******************************************************************/
/* Returns the Legendre polynomial of degree n, evaluated at x.    */
/*******************************************************************/
double legendre(int n, double x) {

    int i;           /* counter */

    double p,        /* Legendre polynomial of order n */
           p_1,      /*    "         "      "    "   n-1*/
           p_2;      /*    "         "      "    "   n-2 */

    p_2 = 1.0;
    p_1 = x;

    if(n >= 2) { 
        for(i=2; i<=n; ++i){
            p = (x*(2.0*i-1.0)*p_1 - (i-1.0)*p_2)/i;
            p_2 = p_1;
            p_1 = p;
        }
        return p;
    } else { 
        if (n == 1) {
            return p_1;
        } else {
            return p_2;
        }
   }
}

/*******************************************************************/
/* Returns the associated Legendre polynomial P_l^m(x).            */
/* Adapted from numerical recipes.                                 */
/*******************************************************************/
double plgndr(int l, int m, double x) {
    double fact, pll, pmm, pmmp1, somx2;
    int i, ll;

    if(m < 0 || m > l || fabs(x) > 1.0) {
        printf("Bad arguments in routine PLGNDR");
    }
    pmm = 1.0;
    if(m > 0) {
        somx2 = sqrt((1.0-x)*(1.0+x));
        fact = 1.0;
        for(i=1; i<=m; ++i) {
            pmm *= -fact*somx2;
            fact += 2.0;
        }
    }
    if(l == m) {
        return pmm;
    } else {
        pmmp1 = x*(2*m+1)*pmm;
        if(l == (m+1)) {
            return pmmp1;
        } else {
            for (ll=(m+2); ll<=l; ll++) {
                pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}

/*******************************************************************/
double rtsec_G(double (*func)(double, double), 
               double Gamma_P, 
               double x1, 
               double x2, 
               double xacc, 
               double ee) {
    int j;
    double fl, f, dx, swap, xl, rts;
 
    fl = (*func)(x1,Gamma_P)-ee;
    f = (*func)(x2,Gamma_P)-ee;

    if(fabs(fl) < fabs(f)) {
        rts = x1;
        xl = x2;
        swap = fl;
        fl = f;
        f = swap;
    } else {
        xl = x1;
        rts = x2;
    }

 
    for(j=1; j<=math::max_iter; ++j) {
        dx = (xl-rts)*f/(f-fl);
        xl = rts;
        fl = f;
        rts += dx;
        f = (*func)(rts,Gamma_P)-ee;

        if(fabs(dx) < xacc || f == 0.0) {
            return rts;
        }
    }
 
    printf("Maximum number of iterations exceeded in rtsec");  
    return 0.0;
}

void print_error(const std::string& error_text) {
    
	std::cerr << "Numerical Recipes run-time error...\n";
	std::cerr << error_text << '\n';
    std::cerr << "...now exiting to system..." << std::endl;
	exit(1);
}
