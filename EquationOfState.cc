#include "EquationOfState.hh"
#include "consts.h"
#include "equil_util.h"
#include <cmath>
#include <cstdio>
#include <stdlib.h>

EquationOfState::EquationOfState() :
    kIsTabulatedEos(false),
    kGammaP(1.0)
{}

EquationOfState::EquationOfState(char eos_file[]) :
    kIsTabulatedEos(true),
    kGammaP(0.0)
{
    loadTabulatedEos(eos_file);    
}

EquationOfState::EquationOfState(double gamma_p) :
    kIsTabulatedEos(false),
    kGammaP(gamma_p)
{}

void EquationOfState::loadTabulatedEos(char eos_file[]) {
    double p,                 /* pressure */
           rho,               /* density */
           h,                 /* enthalpy */
           n0,                /* number density */    
           g;                 /* Gamma */

    FILE *f_eos;              /* pointer to eos_file */
  
    /* OPEN FILE TO READ */
    if((f_eos=fopen(eos_file, "r")) == NULL ) {    
        printf("cannot open file:  %s\n", eos_file); 
        exit(0);
    }

    /* READ NUMBER OF TABULATED POINTS */
    fscanf(f_eos, "%d\n", &_num_tab);

    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS */
    for(int i = 1; i <= _num_tab; ++i) {  
        fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ;
        _log_e_tab[i] = log10(rho*C*C*KSCALE);     /* multiply by C^2 to get */ 
        _log_p_tab[i] = log10(p*KSCALE);           /* energy density. */
        _log_h_tab[i] = log10(h/(C*C));        
        _log_n0_tab[i] = log10(n0);
    }
}
 
double EquationOfState::e_of_rho0(double rho0) const {
    return pow(rho0, kGammaP) / (kGammaP-1.0) + rho0;
}

double EquationOfState::e_of_rho0(double rho0, double Gamma_P) {
    return pow(rho0, Gamma_P) / (Gamma_P-1.0) + rho0;
}

double EquationOfState::e_at_p(double pp, int& n_nearest_pt) const {
    return pow(10.0, interp(_log_p_tab, _log_e_tab, _num_tab, log10(pp), n_nearest_pt));
}

double EquationOfState::e_at_p(double pp) const {
    return pp/(kGammaP-1.0) + pow(pp, 1.0/kGammaP); 
}

double EquationOfState::p_at_e(double ee, int& n_nearest_pt) const {
    return pow(10.0, interp(_log_e_tab, _log_p_tab, _num_tab, log10(ee), n_nearest_pt));
}

double EquationOfState::p_at_h(double hh, int& n_nearest_pt) const {
    return pow(10.0, interp(_log_h_tab, _log_p_tab, _num_tab, log10(hh), n_nearest_pt));
}

double EquationOfState::h_at_p(double pp, int& n_nearest_pt) const {
    return pow(10.0, interp(_log_p_tab, _log_h_tab, _num_tab, log10(pp), n_nearest_pt));
}

double EquationOfState::n0_at_e(double ee, int& n_nearest_pt) const {
    return pow(10.0, interp(_log_e_tab, _log_n0_tab, _num_tab, log10(ee), n_nearest_pt));
}

bool EquationOfState::isTabulatedEos() const {
    return kIsTabulatedEos;
}

double EquationOfState::getGammaP() const {
    return kGammaP;
}

int EquationOfState::getNumTab() const {
    return _num_tab;
}
