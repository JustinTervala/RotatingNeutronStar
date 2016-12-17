#include "RotatingNeutronStar.hh"


RotatingNeutonStar::RotatingNeutronStar(const char eos_file[], double central_energy_density) :
        eos(eos_file),
        e_center(central_energy_density),
        e_surface(7.8*C*C*KSCALE),
        p_surface(1.01e8*KSCALE),
        enthalpy_min(1.0/(C*C))   
{
    initalizeRns();
}

RotatingNeutonStar::RotatingNeutronStar(double gamma_p, double central_energy_density) :
        eos(gamma_p),
        e_center(central_energy_density),
        e_surface(0.0),
        p_surface(0.0),
        enthalpy_min(0.0)   
{
    initializeRns();
}

void RotatingNeutronStar::initializeRns() {
    make_grid();
    make_center();
    sphere();
}

void RotatingNeutronStar::make_grid() {
    int m, s;                         /* counters */
    
    for(s = 1; s <= SDIV; ++s) {
        s_gp[s] = SMAX*(s-1.0)/(SDIV-1.0);
    }
    /* s_gp[1] = 0.0     corresponds to the center of the star
       s_gp[SDIV] = SMAX corresponds to infinity */

    for(m = 1; m<= MDIV; ++m) { 
         mu[m] = (m-1.0)/(MDIV-1.0);
    }
    /* mu[1] = 0.0    corresponds to the plane of the equator 
       mu[MDIV] = 1.0 corresponds to the axis of symmetry */

    /* s_gp[0] and mu[0] are not used by the program */
}

//TODO: Thsi could be done in initializer list?
void RotatingNeutronStar::make_center() {

    int n_nearest = eos.getNumTab()/2; 

    if(eos.isTabulatedEos()) {
        p_center = eos.p_at_e(e_center, n_nearest);
        h_center = eos.h_at_p(p_center, n_nearest);
    } else {
        double rho0_center = rtsec_G(EquationOfState::e_of_rho0, eos.getGammaP(), 0.0, e_center, DBL_EPSILON, e_center );
        p_center = pow(rho0_center, eos.getGammaP());
        h_center = log((e_center+p_center)/rho0_center);
    }
}

void RotatingNeutronStar::sphere() {
    int s, m, n_nearest;

    double r_is_s,
           r_is_final,
           r_final, 
           m_final,
           lambda_s,
           nu_s,
           r_is_gp[RDIV+1],
           lambda_gp[RDIV+1],
           nu_gp[RDIV+1],
           gama_mu_0[SDIV+1],
           rho_mu_0[SDIV+1],
           gama_eq,
           rho_eq,
           s_e=0.5;

    /* The function TOV integrates the TOV equations. The function
       can be found in the file equil.c */

    TOV(1, e_center, p_center, p_surface, e_surface, eos, 
        r_is_gp, lambda_gp, nu_gp, r_is_final, r_final, m_final);

    TOV(2, e_center, p_center, p_surface, e_surface, eos,
        r_is_gp, lambda_gp, nu_gp, r_is_final, r_final, m_final);

    TOV(3, e_center, p_center, p_surface, e_surface, eos, 
        r_is_gp, lambda_gp, nu_gp, r_is_final, r_final, m_final);

    n_nearest = RDIV/2;
    for(s=1; s<=SDIV; ++s) {
        r_is_s = r_is_final*(s_gp[s]/(1.0-s_gp[s]));

        if(r_is_s < r_is_final) {
            lambda_s = interp(r_is_gp, lambda_gp, RDIV, r_is_s, n_nearest);
            nu_s = interp(r_is_gp, nu_gp, RDIV, r_is_s, n_nearest);
        } else {
            lambda_s = 2.0*log(1.0+m_final/(2.0*r_is_s));
            nu_s = log((1.0-m_final/(2.0*r_is_s))/(1.0+m_final/(2*r_is_s)));
        }

        metric[s][1].gama = nu_s+lambda_s;
        metric[s][1].rho = nu_s-lambda_s;

        for(m=1; m<=MDIV; ++m) {
            metric[s][m].rho = metric[s][1].rho;
            metric[s][m].gama = metric[s][1].gama;        
            metric[s][m].omega = 0.0; 
            metric[s][m].alpha = (metric[s][1].gama-metric[s][1].rho)/2.0;
        }
 
        rho_mu_0[s] = metric[s][1].rho;                     /* rho at \mu=0 */
        gama_mu_0[s] = metric[s][1].gama;                   /* gama at \mu=0 */

    }

    n_nearest = SDIV/2;
    gama_eq = interp(s_gp, gama_mu_0, SDIV, s_e, n_nearest); /* gama at equator */
    rho_eq = interp(s_gp, rho_mu_0, SDIV, s_e, n_nearest);   /* rho at equator */
 
    r_e = r_final*exp(0.5*(rho_eq-gama_eq)); 
}



void RotatingNeutronStar::TOV(int i_check, 
                              double r_is_gp[RDIV+1], 
                              double lambda_gp[RDIV+1], 
                              double nu_gp[RDIV+1], 
                              double &r_is_final, 
                              double &r_final, 
                              double &m_final) {

    int i = 2, n_nearest;

    double r,                           /* radius */
           r_is,                        /* isotropic radial coordinate */
           r_is_est,                    /* estimate on final isotr. radius */ 
           r_is_check,                  /*                      */    
           dr_is_save,                  /* r_is saving interval */  
           rho0,
           e_d,                         /* density */
           p,                           /* pressure */
           h,                           /* stepsize during integration */
           m,                           /* mass   */
           nu_s,
           hh,
           a1,a2,a3,a4,b1,b2,b3,b4,     /* coeff. in Runge-Kutta equations */
           c1,c2,c3,c4,
           k_rescale, 
           r_gp[RDIV+1],
           m_gp[RDIV+1],
           e_d_gp[RDIV+1];   

    if(i_check == 1) {
        if(eos.isTabulatedEos()) {
            r_is_est=1.5e6/sqrt(KAPPA);
        } else {
            r_is_est=2.0*sqrt(eos.getGammaP()/(4.0*PI*(eos.getGammaP()-1.0)))*pow(e_center, (eos.getGammaP()-2.0)/2.0);
        }
        h=r_is_est/100;     
    } else {
        r_is_est= r_is_final;
        h = r_is_est/10000;   
        dr_is_save = r_is_final/RDIV;
        r_is_check = dr_is_save;
    }

    r_is = 0.0;                            /* initial isotropic radius */
    r = 0.0;                               /* initial radius */
    m = 0.0;                               /* initial mass */ 
    p = p_center;                          /* initial pressure */ 

    r_is_gp[1] = 0.0;
    r_gp[1] = 0.0;
    m_gp[1] = 0.0;
    lambda_gp[1] = 0.0;
    e_d_gp[1] = e_center; 

    n_nearest = eos.getNumTab()/2;

    while(p >= p_surface) { 
 
        e_d = eos.e_at_p(p, n_nearest);

        if((i_check == 3) && (r_is > r_is_check) && (i <= RDIV)) {
            r_is_gp[i] = r_is;
            r_gp[i] = r;
            m_gp[i] = m;
            e_d_gp[i] = e_d; 
            ++i;   
            r_is_check += dr_is_save;
        }    
       
        r_is_final = r_is;
        r_final = r;
        m_final = m;

        a1 = dr_dr_is(r_is, r, m);

        //TODO: Look at how n_nearest changes. May be best to group similar n_nearest together
        b1 = dm_dr_is(r_is, r, m, p, e_center, p_surface, eos, n_nearest);
        c1 = dp_dr_is(r_is,r,m,p, e_center, p_surface, eos, n_nearest);

        a2 = dr_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0);

        b2 = dm_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0, e_center, 
                      p_surface, eos, n_nearest);

        c2 = dp_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0, e_center, 
                      p_surface, eos, n_nearest);  

        a3 = dr_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0);

        b3 = dm_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0, e_center, 
                      p_surface, eos, n_nearest);

        c3 = dp_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0, e_center, 
                      p_surface, eos, n_nearest); 

        a4 = dr_dr_is(r_is+h, r+h*a3, m+h*b3);

        b4 = dm_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3, e_center, p_surface, eos, n_nearest); 

        c4 = dp_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3, e_center, p_surface, eos, n_nearest); 

        r += (h/6.0)*(a1+2*a2+2*a3+a4);
        m += (h/6.0)*(b1+2*b2+2*b3+b4);
        p += (h/6.0)*(c1+2*c2+2*c3+c4);

        r_is += h;

    }

    r_is_gp[RDIV] = r_is_final;
    r_gp[RDIV] = r_final;
    m_gp[RDIV] = m_final;

    /* Rescale r_is and compute lambda */

    if(i_check == 3) {
        k_rescale = 0.5*(r_final/r_is_final)*(1.0-m_final/r_final+
                    sqrt(1.0-2.0*m_final/r_final));
 
        r_is_final *= k_rescale;
 
        nu_s = log((1.0-m_final/(2.0*r_is_final))/(1.0+m_final/
               (2.0*r_is_final)));

        for(i=1; i<=RDIV; ++i) {
            r_is_gp[i] *= k_rescale;
 
            if(i == 1) {
                lambda_gp[1] = log(1.0/k_rescale);
            } else {
                lambda_gp[i] = log(r_gp[i]/r_is_gp[i]); 
            }
            if(e_d_gp[i] < e_surface) {
                hh=0.0;
            } else { 
                if(eos.isTabulatedEos()) {
                    p = eos.p_at_e(e_d_gp[i], n_nearest);
                    hh = eos.h_at_p(p, n_nearest);
                } else { 
                    rho0 = rtsec_G(EquationOfState::e_of_rho0, eos.getGammaP(), 0.0, e_d_gp[i], DBL_EPSILON, e_d_gp[i]);
                    p = pow(rho0, eos.getGammaP());
                    hh = log((e_d_gp[i]+p)/rho0);
                }
            }
 
            nu_gp[i] = nu_s-hh;
        }
        nu_gp[RDIV]=nu_s;
    }
}


double RotatingNeutronStar::dm_dr_is(double r_is, 
                                     double r, 
                                     double m, 
                                     double p, 
                                     int &n_nearest_pt) {

    double dmdr, e_d;

    if(p < p_surface) { 
        e_d = 0.0;
    } else { 
        e_d = eos.e_at_p(p, n_nearest_pt);
    } 
    if(r_is < RMIN) { 
        dmdr=4.0*PI*e_center*r*r*(1.0+4.0*PI*e_center*r*r/3.0);
    } else {
        dmdr=4.0*PI*e_d*r*r*r*sqrt(1.0-2.0*m/r)/r_is;
    }
    return dmdr;
}


double RotatingNeutronStar::dp_dr_is(double r_is, 
                                     double r, 
                                     double m, 
                                     double p,
                                     int &n_nearest_pt) {

    double dpdr, e_d; 

    if(p < p_surface) { 
        e_d=0.0;
    } else {        
        e_d=eos.e_at_p(p, n_nearest_pt);
    }
    if(r_is < RMIN) {
        dpdr = -4.0*PI*(e_center+p)*(e_center+3.0*p)*r*(1.0+4.0*e_center*r*r/3.0)/3.0;
    } else {
        dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*p)/(r*r_is*sqrt(1.0-2.0*m/r));
    }
    return dpdr;
}

double RotatingneutronStar::dr_dr_is(double r_is, double r, double m) {
    if(r_is < RMIN) {
        drdris=1.0;
    } else {
        drdris=(r/r_is)*sqrt(1.0-2.0*m/r);
    }
    return drdris;
}


void RotatingNeutronStar::mass_radius(double r_ratio) {
    int s,
        m,
        n_nearest;

 
    double //**rho_0, /*rest mass density*/
           //**velocity,
           gama_equator,              /* gama at equator */
           rho_equator,               /* rho at equator */
           omega_equator,             /* omega at equator */
           s1,
           s_1,
           d_gama_s,
           d_rho_s,
           d_omega_s,
           sqrt_v,
           D_m[SDIV+1],               /* int. quantity for M */
           D_m_0[SDIV+1],             /* int. quantity for M_0 */ 
           D_J[SDIV+1],               /* int. quantity for J */
           s_e,                 
           d_o_e[SDIV+1],
           d_g_e[SDIV+1],
           d_r_e[SDIV+1],
           d_v_e[SDIV+1],
           doe,
           dge, 
           dre,
           dve,
           vek,     
           gama_mu_0[SDIV+1],                   
           rho_mu_0[SDIV+1],                    
           omega_mu_0[SDIV+1],
           J,
           r_p,
           s_p;                 
        
    r_p = r_ratio*r_e;                              /* radius at pole */
    s_p = r_p/(r_p+r_e);                            /* s-coordinate at pole */
    s_e = 0.5;
    
    native_matrix<double, SDIV+1, MDIV+1> rho_0;
    matrix<double, SDIV+1, MDIV+1> velocity = {{0.0}};
    for(s = 1; s <= SDIV; ++s) {               
        gama_mu_0[s]=metric[s][1].gama;                   
        rho_mu_0[s]=metric[s][1].rho;                                                    
    }

    n_nearest = SDIV/2;
    gama_equator = interp(s_gp,gama_mu_0,SDIV,s_e, n_nearest);  
    rho_equator = interp(s_gp,rho_mu_0,SDIV,s_e, n_nearest);   

    /* Circumferential radius */

    if(eos.isTabulatedEos()) {
        R_e = sqrt(KAPPA)*r_e*exp((gama_equator-rho_equator)/2.0);
    } else {
        R_e = r_e*exp((gama_equator-rho_equator)/2.0);
    }
    
    /* Masses and angular momentum */
 
    Mass = 0.0;              /* initialize */
    Mass_0 = 0.0;
    J = 0.0;

   /* CALCULATE THE REST MASS DENSITY */
    if(eos.isTabulatedEos()) {
        n_nearest = eos.getNumTab()/2;
        for(s=1; s<=SDIV; ++s) {
            for(m=1; m<=MDIV; ++m) {
                if(energy[s][m] > e_surface) {
                    rho_0[s][m] = eos.n0_at_e(energy[s][m], n_nearest)*MB*KSCALE*SQ(C);
                } else {
                    rho_0[s][m] = 0.0;
                }
            }
        }       
    } else {
        for(s=1; s<=SDIV; ++s) {
            for(m=1; m<=MDIV; ++m) {
                rho_0[s][m] = (energy[s][m]+pressure[s][m])*exp(-enthalpy[s][m]);
            }
        }
    }

    for(s=1; s<=SDIV; ++s) {
        D_m[s] = 0.0;           /* initialize */
        D_m_0[s] = 0.0;
        D_J[s] = 0.0;

        for(m=1; m<=MDIV-2; m+=2) {
            D_m[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*metric[s][m].alpha+metric[s][m].gama)*
                      (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
                      (1.0+velocity_sq[s][m]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
                      (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_e*metric[s][m].omega*
                      exp(-metric[s][m].rho)) + 2.0*pressure[s][m])

                    + 4.0*exp(2.0*metric[s][m+1].alpha+metric[s][m+1].gama)*
                      (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
                      (1.0+velocity_sq[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
                      (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_e*metric[s][m+1].omega*
                      exp(-metric[s][m+1].rho)) + 2.0*pressure[s][m+1]) 

                    + exp(2.0*metric[s][m+2].alpha+metric[s][m+2].gama)*
                      (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
                      (1.0+velocity_sq[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
                      (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_e*metric[s][m+2].omega*
                      exp(-metric[s][m+2].rho)) + 2.0*pressure[s][m+2]));    

            D_m_0[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*metric[s][m].alpha+(metric[s][m].gama
                        -metric[s][m].rho)/2.0)*rho_0[s][m]/sqrt(1.0-velocity_sq[s][m])

                     + 4.0* exp(2.0*metric[s][m+1].alpha+(metric[s][m+1].gama
                       -metric[s][m+1].rho)/2.0)*rho_0[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])
         
                     + exp(2.0*metric[s][m+2].alpha+(metric[s][m+2].gama
                       -metric[s][m+2].rho)/2.0)*rho_0[s][m+2]/sqrt(1.0-velocity_sq[s][m+2]));
  
            D_J[s] += (1.0/(3.0*(MDIV-1)))*( sqrt(1.0-mu[m]*mu[m])*
                      exp(2.0*metric[s][m].alpha+metric[s][m].gama-metric[s][m].rho)*(energy[s][m]
                      +pressure[s][m])*sqrt(velocity_sq[s][m])/(1.0-velocity_sq[s][m])
  
                   + 4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
                     exp(2.0*metric[s][m+1].alpha+metric[s][m+1].gama-metric[s][m+1].rho)*(energy[s][m+1]
                     +pressure[s][m+1])*sqrt(velocity_sq[s][m+1])/
                     (1.0-velocity_sq[s][m+1])

                   + sqrt(1.0-mu[m+2]*mu[m+2])*
                     exp(2.0*metric[s][m+2].alpha+metric[s][m+2].gama-metric[s][m+2].rho)*(energy[s][m+2]
                     +pressure[s][m+2])*sqrt(velocity_sq[s][m+2])/
                     (1.0-velocity_sq[s][m+2]));
        }
    }

    for(s=1; s<=SDIV-2; s+=2) { 
        Mass += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]), 4.0)*
                   D_m[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]), 4.0)*D_m[s+1]
                   +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]), 4.0)*D_m[s+2]);

        Mass_0 += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]), 4.0)*
                     D_m_0[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]), 4.0)*D_m_0[s+1]
                     +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]), 4.0)*D_m_0[s+2]);
 
        J += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s], 3.0)/pow(1.0-s_gp[s], 5.0))*
             D_J[s]+ 4.0*(pow(s_gp[s+1], 3.0)/pow(1.0-s_gp[s+1], 5.0))*
             D_J[s+1] + (pow(s_gp[s+2], 3.0)/pow(1.0-s_gp[s+2], 5.0))*
             D_J[s+2]);

    }
   
    if(eos.isTabulatedEos()) {
        Mass *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e, 3.0)/G;
        Mass_0 *= 4*PI*sqrt(KAPPA)*C*C*pow(r_e, 3.0)/G;
    } else {
        Mass *= 4*PI*pow(r_e, 3.0);
        Mass_0 *= 4*PI*pow(r_e, 3.0);
    }
 
    if(r_ratio == 1.0) { 
        J = 0.0; 
    } else {    
        if(eos.isTabulatedEos()) {
            J *= 4*PI*KAPPA*C*C*C*pow(r_e, 4.0)/G;
        } else { 
            J *= 4*PI*pow(r_e,4.0);
        }
    }

    ang_mom = J;


  /* Compute the velocities of co-rotating and counter-rotating particles
    with respect to a ZAMO     */

    for(s=1+(SDIV-1)/2; s<=SDIV; ++s) {
        s1 = s_gp[s]*(1.0-s_gp[s]);
        s_1 = 1.0-s_gp[s];
        
        d_rho_s = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::rho, s, 1);
        d_gama_s = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::gama, s, 1);
        d_omega_s = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::omega, s, 1);

        sqrt_v = exp(-2.0*metric[s][1].rho)*r_e*r_e*pow(s_gp[s], 4.0)*pow(d_omega_s,2.0) 
                 + 2*s1*(d_gama_s+d_rho_s)+s1*s1*(d_gama_s*d_gama_s-d_rho_s*d_rho_s);

        if(sqrt_v > 0.0) {
            sqrt_v = sqrt(sqrt_v);
        } else {
            sqrt_v=0.0;
        }

        v_plus[s]=(exp(-metric[s][1].rho)*r_e*s_gp[s]*s_gp[s]*d_omega_s + sqrt_v)/
                  (2.0+s1*(d_gama_s-d_rho_s));

        v_minus[s]=(exp(-metric[s][1].rho)*r_e*s_gp[s]*s_gp[s]*d_omega_s - sqrt_v)/
                   (2.0+s1*(d_gama_s-d_rho_s));
    }


    /* Kepler angular velocity */

    for(s=1; s<=SDIV; ++s) { 
        d_r_e[s] = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::rho, s, 1);
        d_g_e[s] = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::gama, s, 1);
        d_o_e[s] = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::omega, s, 1);
        d_v_e[s] = deriv_s<double, SDIV+1, MDIV+1>(velocity, s, 1);
        /* Value of omega on the equatorial plane*/
        omega_mu_0[s] = metric[s][1].omega;
    }

    n_nearest=SDIV/2; 
    doe=interp(s_gp, d_o_e, SDIV, s_e, n_nearest);
    dge=interp(s_gp, d_g_e, SDIV, s_e, n_nearest);
    dre=interp(s_gp, d_r_e, SDIV, s_e, n_nearest);
    dve=interp(s_gp, d_v_e, SDIV, s_e, n_nearest);

    vek = (doe/(8.0+dge-dre))*r_e*exp(-rho_equator) + sqrt(((dge+dre)/(8.0+dge
          -dre)) + pow((doe/(8.0+dge-dre))*r_e*exp(-rho_equator),2.0));


    if (r_ratio ==1.0) {
        omega_equator = 0.0;
    } else {
        omega_equator = interp(s_gp, omega_mu_0, SDIV, s_e, n_nearest);
    }


    if(eos.isTabulatedEos()) {
        Omega_K = (C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_e);
    } else { 
        Omega_K = omega_equator + vek*exp(rho_equator)/r_e;
    }
}


void RotatingNeutronStar::spin(double r_ratio) {
    int m,                      /* counter */
        s,                      /* counter */
        n,                      /* counter */
        k,                      /* counter */
        n_of_it = 0,            /* number of iterations */
        n_nearest,
        print_dif = 0,
        i,
        j;
    int a_check = 0;
    double sum_rho = 0.0,         /* intermediate sum in eqn for rho */
           sum_gama = 0.0,        /* intermediate sum in eqn for gama */
           sum_omega=0.0,       /* intermediate sum in eqn for omega */
           r_e_old,             /* equatorial radius in previus cycle */
           dif = 1.0,             /* difference | r_e_old - r_e | */
           d_gama_s,            /* derivative of gama w.r.t. s */
           d_gama_m,            /* derivative of gama w.r.t. m */
           d_rho_s,             /* derivative of rho w.r.t. s */
           d_rho_m,             /* derivative of rho w.r.t. m */
           d_omega_s,           /* derivative of omega w.r.t. s */
           d_omega_m,           /* derivative of omega w.r.t. m */
           d_gama_ss,           /* 2nd derivative of gama w.r.t. s */
           d_gama_mm,           /* 2nd derivative of gama w.r.t. m */
           d_gama_sm,           /* derivative of gama w.r.t. m and s */
           temp1,                /* temporary term in da_dm */ 
           temp2, 
           temp3,
           temp4,
           temp5,
           temp6,
           temp7,
           temp8,
           m1,                  
           s1,
           s2,
           ea,
           rsm,
           gsm,
           omsm,
           esm,
           psm,
           v2sm,
           mum,
           sgp,
           s_1,
           e_gsm,
           e_rsm, 
           rho0sm,
           term_in_Omega_h,
           r_p,
           s_p,
           gama_pole_h,                  /* gama^hat at pole */  
           gama_center_h,                /* gama^hat at center */
           gama_equator_h,               /* gama^hat at equator */
           rho_pole_h,                   /* rho^hat at pole */ 
           rho_center_h,                 /* rho^hat at center */
           rho_equator_h,                /* rho^hat at equator */ 
           omega_equator_h,              /* omega^hat at equator */         
           gama_mu_1[SDIV+1],            /* gama at \mu=1 */
           gama_mu_0[SDIV+1],            /* gama at \mu=0 */
           rho_mu_1[SDIV+1],             /* rho at \mu=1 */
           rho_mu_0[SDIV+1],             /* rho at \mu=0 */
           omega_mu_0[SDIV+1],           /* omega at \mu=0 */
           s_e = 0.5,
           Omega_h,
           sk,
           sj,
           sk1,
           sj1,
           r_e;

  
    r_e = r_e_new;
    native_matrix<RhoGamaOmega, SDIV+1, MDIV+1> S_metric;
    native_matrix<RhoGamaOmega, LMAX+2, SDIV+1> D1_metric;
    native_matrix<RhoGamaOmega, SDIV+1, LMAX+2> D2_metric;
    matrix<double, SDIV+1, MDIV+1> da_dm = {{0.0}};
    matrix<double, SDIV+1, MDIV+1> dgds = {{0.0}};
    matrix<double, SDIV+1, MDIV+1> dgdm = {{0.0}};

    while(dif > control.accuracy || n_of_it < 2) { 

        if(print_dif != 0) {
            printf("%4.3e\n", dif);
        }
 
        /* Rescale potentials and construct arrays with the potentials along
        | the equatorial and polar directions.
        */        

        for(s=1; s<=SDIV; ++s) {
            for(m=1; m<=MDIV; ++m) {
                metric[s][m].rho /= SQ(r_e);
                metric[s][m].gama /= SQ(r_e); 
                metric[s][m].alpha /= SQ(r_e);
                metric[s][m].omega *= r_e;
            }
            rho_mu_0[s] = metric[s][1].rho;     
            gama_mu_0[s] = metric[s][1].gama;   
            omega_mu_0[s] = metric[s][1].omega; 
            rho_mu_1[s] = metric[s][MDIV].rho;  
            gama_mu_1[s] = metric[s][MDIV].gama;
        }
 
        /* Compute new r_e. */ 

        r_e_old = r_e;
        r_p = r_ratio*r_e;                          
        s_p = r_p/(r_p+r_e);                        
  
        n_nearest = SDIV/2;
        gama_pole_h = interp(s_gp, gama_mu_1, SDIV, s_p, n_nearest); 
        gama_equator_h = interp(s_gp, gama_mu_0, SDIV, s_e, n_nearest);
        gama_center_h = metric[1][1].gama;                    
  
        rho_pole_h = interp(s_gp, rho_mu_1, SDIV, s_p, n_nearest);   
        rho_equator_h = interp(s_gp, rho_mu_0, SDIV, s_e, n_nearest);
        rho_center_h = metric[1][1].rho;                      
 
        r_e = sqrt(2*h_center/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h));

        /* Compute angular velocity Omega. */
 
        if(r_ratio == 1.0) {
            Omega_h = 0.0;
            omega_equator_h = 0.0;
        } else {
            omega_equator_h = interp(s_gp, omega_mu_0, SDIV, s_e, n_nearest);
            term_in_Omega_h = 1.0-exp(SQ(r_e)*(gama_pole_h+rho_pole_h
                                              -gama_equator_h-rho_equator_h));
            if(term_in_Omega_h >= 0.0) { 
               Omega_h = omega_equator_h + exp(SQ(r_e)*rho_equator_h)
                                            *sqrt(term_in_Omega_h);
            } else {
               Omega_h = 0.0;
            }
        }
 
        /* Compute velocity, energy density and pressure. */
        struct timespec vep_start, vep_stop;
        clock_gettime(CLOCK_MONOTONIC, &vep_start);
 
        n_nearest = eos.getNumTab()/2; 

        for(s=1; s<=SDIV; ++s) {
            sgp = s_gp[s];

            for(m=1; m<=MDIV; ++m) {
                rsm = metric[s][m].rho;
            
                if(r_ratio != 1.0) { 
                    velocity_sq[s][m] = SQ((Omega_h-metric[s][m].omega)*(sgp/(1.0-sgp))
                                           *trig.sin_theta[m]*exp(-rsm*SQ(r_e)));
                } else {
                    velocity_sq[s][m] = 0.0;
		}
                if(velocity_sq[s][m] >= 1.0) { 
                    velocity_sq[s][m] = 0.0;
                }
                enthalpy[s][m] = enthalpy_min + 0.5*(SQ(r_e)*(gama_pole_h+rho_pole_h
                                 -metric[s][m].gama-rsm)-log(1.0-velocity_sq[s][m]));
  
                if((enthalpy[s][m] <= enthalpy_min) || (sgp > s_e)) {
                    pressure[s][m] = 0.0;
                    energy[s][m] = 0.0; 
                } else { 
                    if(eos.isTabulatedEos()) {
                        pressure[s][m] = eos.p_at_h(enthalpy[s][m], n_nearest);
                        energy[s][m] = eos.e_at_p(pressure[s][m], n_nearest);
                    } else {
                        rho0sm = pow(((eos.getGammaP()-1.0)/eos.getGammaP())
                                     *(exp(enthalpy[s][m])-1.0), 1.0/(eos.getGammaP()-1.0));
 
                        pressure[s][m] = pow(rho0sm, eos.getGammaP());

                        energy[s][m] = pressure[s][m]/(eos.getGammaP()-1.0)+rho0sm;
                    }
                }  

                /* Rescale back metric potentials (except omega) */

                metric[s][m].rho *= SQ(r_e);
                metric[s][m].gama *= SQ(r_e);
                metric[s][m].alpha *= SQ(r_e);
            }
        }
        clock_gettime(CLOCK_MONOTONIC, &vep_stop);
        printf("spin(), vep: %ld\n", getElapsedTimeNs(vep_start, vep_stop));

        /* Compute metric potentials */
        struct timespec metric_start, metric_stop;
        clock_gettime(CLOCK_MONOTONIC, &metric_start);

        for(s=1; s<=SDIV; ++s) {
            for(m=1; m<=MDIV; ++m) {
                rsm = metric[s][m].rho;
                gsm = metric[s][m].gama;
                omsm = metric[s][m].omega;
                esm = energy[s][m];
                psm = pressure[s][m];
                e_gsm = exp(0.5*gsm);
                e_rsm = exp(-rsm);
                v2sm = velocity_sq[s][m];
                mum = mu[m];            
                m1 = 1.0-SQ(mum);
                sgp = s_gp[s];
                s_1 = 1.0-sgp;
                s1 = sgp*s_1;
                s2 = SQ(sgp/s_1);  

                ea = 16.0*PI*exp(2.0*metric[s][m].alpha)*SQ(r_e);
 
                if(s == 1) {
                    d_gama_s = 0.0;
                    d_gama_m = 0.0;
                    d_rho_s = 0.0;
                    d_rho_m = 0.0;
                    d_omega_s = 0.0;
                    d_omega_m = 0.0;
                } else {
                    d_rho_s = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::rho, s, m);
                    d_rho_m = deriv_m<SDIV+1, MDIV+1>(metric, &Metric::rho, s, m);
                    d_gama_s = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::gama, s, m);
                    d_gama_m = deriv_m<SDIV+1, MDIV+1>(metric, &Metric::gama, s, m);
                    d_omega_s = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::omega, s, m);
                    d_omega_m = deriv_m<SDIV+1, MDIV+1>(metric, &Metric::omega, s, m);
                }      

                S_metric[s][m].rho = e_gsm*(0.5*ea*(esm + psm)*s2*(1.0+v2sm)/(1.0-v2sm)
                                + s2*m1*SQ(e_rsm)*(SQ(s1*d_omega_s) 
                                + m1*SQ(d_omega_m))
                                + s1*d_gama_s - mum*d_gama_m + 0.5*rsm*(ea*psm*s2  
                                - s1*d_gama_s*(0.5*s1*d_gama_s+1.0) 
                                - d_gama_m*(0.5*m1*d_gama_m-mum)));

                S_metric[s][m].gama = e_gsm*(ea*psm*s2 + 0.5*gsm*(ea*psm*s2 - 0.5*SQ(s1
                                *d_gama_s) - 0.5*m1*SQ(d_gama_m)));

                S_metric[s][m].omega = e_gsm*e_rsm*( -ea*(Omega_h-omsm)*(esm+psm)
                                *s2/(1.0-v2sm) + omsm*( -0.5*ea*(((1.0+v2sm)*esm 
                                + 2.0*v2sm*psm)/(1.0-v2sm))*s2 
                                - s1*(2*d_rho_s+0.5*d_gama_s)
                                + mum*(2*d_rho_m+0.5*d_gama_m) + 0.25*SQ(s1)*(4
                                *SQ(d_rho_s)-SQ(d_gama_s)) + 0.25*m1*(4*SQ(d_rho_m)
                                - SQ(d_gama_m)) - m1*SQ(e_rsm)*(SQ(SQ(sgp)*d_omega_s)
                                + s2*m1*SQ(d_omega_m))));
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &metric_stop);
        printf("spin(), metric: %ld\n", getElapsedTimeNs(metric_start, metric_stop));
        /* ANGULAR INTEGRATION */
        struct timespec ang_start, ang_stop;
        clock_gettime(CLOCK_MONOTONIC, &ang_start);
   
        n = 0;
        for(k=1; k<=SDIV; ++k) {      
            for(m=1; m<=MDIV-2; m+=2) {
                sum_rho += trig.P_2n[m][n+1]*S_metric[k][m].rho
                           + 4.0*trig.P_2n[m+1][n+1]*S_metric[k][m+1].rho 
                           + trig.P_2n[m+2][n+1]*S_metric[k][m+2].rho;
            }

            D1_metric[n+1][k].rho = (DM/3.0)*sum_rho;
            D1_metric[n+1][k].gama = 0.0;
            D1_metric[n+1][k].omega = 0.0;
            sum_rho = 0.0;
        }

        for(n=1; n<=LMAX; ++n) {
            for(k=1; k<=SDIV; ++k) {      
                for(m=1; m<=MDIV-2; m+=2) {

                    sum_rho += trig.P_2n_t[n+1][m]*S_metric[k][m].rho
                               + 4.0*trig.P_2n_t[n+1][m+1]*S_metric[k][m+1].rho 
                               + trig.P_2n_t[n+1][m+2]*S_metric[k][m+2].rho;
                       
                    sum_gama += trig.sin_2n_1_theta_t[n][m]*S_metric[k][m].gama
                                +4.0*trig.sin_2n_1_theta_t[n][m+1]*S_metric[k][m+1].gama
                                +trig.sin_2n_1_theta_t[n][m+2]*S_metric[k][m+2].gama;
  
                    sum_omega += trig.sin_theta[m]*trig.P1_2n_1_t[n+1][m]*S_metric[k][m].omega
                                 +4.0*trig.sin_theta[m+1]*trig.P1_2n_1_t[n+1][m+1]*S_metric[k][m+1].omega
                                 +trig.sin_theta[m+2]*trig.P1_2n_1_t[n+1][m+2]*S_metric[k][m+2].omega;
                }
                D1_metric[n+1][k].rho = (DM/3.0)*sum_rho;
                D1_metric[n+1][k].gama = (DM/3.0)*sum_gama;
                D1_metric[n+1][k].omega = (DM/3.0)*sum_omega;
                sum_rho = 0.0;
                sum_gama = 0.0;
                sum_omega = 0.0;
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &ang_stop);
        printf("spin(), ang: %ld\n", getElapsedTimeNs(ang_start, ang_stop));
        /* RADIAL INTEGRATION */
        struct timespec rad_start, rad_stop;
        clock_gettime(CLOCK_MONOTONIC, &rad_start);



        n = 0;
        for(s=1; s<=SDIV; ++s) {
            for(k=1; k<=SDIV-2; k+=2) { 
                sum_rho += ( trig.f_rho[s][n+1][k]*D1_metric[n+1][k].rho 
                           + 4.0*trig.f_rho[s][n+1][k+1]*D1_metric[n+1][k+1].rho
                            + trig.f_rho[s][n+1][k+2]*D1_metric[n+1][k+2].rho);
            }
            D2_metric[s][n+1].rho = (DS/3.0)*sum_rho;
            D2_metric[s][n+1].gama = 0.0;
            D2_metric[s][n+1].omega = 0.0;
            sum_rho = 0.0;
        }  
 
        for(s=1; s<=SDIV; ++s) {
            for(n=1; n<=LMAX; ++n) {
                for(k=1; k<=SDIV-2; k+=2) { 
                    sum_rho += trig.f_rho[s][n+1][k]*D1_metric[n+1][k].rho 
                               + 4.0*trig.f_rho[s][n+1][k+1]*D1_metric[n+1][k+1].rho
                               + trig.f_rho[s][n+1][k+2]*D1_metric[n+1][k+2].rho;
 
                    sum_gama += trig.f_gama[s][n+1][k]*D1_metric[n+1][k].gama 
                                + 4.0*trig.f_gama[s][n+1][k+1]*D1_metric[n+1][k+1].gama
                                + trig.f_gama[s][n+1][k+2]*D1_metric[n+1][k+2].gama;
     
                    if(k < s && k+2 <= s) {
                        sum_omega += trig.f_rho[s][n+1][k]*D1_metric[n+1][k].omega 
                                     + 4.0*trig.f_rho[s][n+1][k+1]*D1_metric[n+1][k+1].omega
                                     + trig.f_rho[s][n+1][k+2]*D1_metric[n+1][k+2].omega;
                    } else {
                        if(k >= s) {
                            sum_omega += trig.f_gama[s][n+1][k]*D1_metric[n+1][k].omega 
                                         + 4.0*trig.f_gama[s][n+1][k+1]*D1_metric[n+1][k+1].omega
                                         + trig.f_gama[s][n+1][k+2]*D1_metric[n+1][k+2].omega;
                        } else {
                            sum_omega += trig.f_rho[s][n+1][k]*D1_metric[n+1][k].omega 
                                         + 4.0*trig.f_rho[s][n+1][k+1]*D1_metric[n+1][k+1].omega
                                         + trig.f_gama[s][n+1][k+2]*D1_metric[n+1][k+2].omega;
                        }
                    }
                }
                D2_metric[s][n+1].rho = (DS/3.0)*sum_rho;
                D2_metric[s][n+1].gama = (DS/3.0)*sum_gama;
                D2_metric[s][n+1].omega = (DS/3.0)*sum_omega;
                sum_rho = 0.0;
                sum_gama = 0.0;
                sum_omega = 0.0;
            }
        }   
 
        clock_gettime(CLOCK_MONOTONIC, &rad_stop);
        printf("spin(), rad: %ld\n", getElapsedTimeNs(rad_start, rad_stop));

        /* SUMMATION OF COEFFICIENTS */
        struct timespec coeff_start, coeff_stop;
        clock_gettime(CLOCK_MONOTONIC, &coeff_start);

        for(s=1; s<=SDIV; ++s) {
            for(m=1; m<=MDIV; ++m) {

                rsm = metric[s][m].rho;
                gsm = metric[s][m].gama;
                omsm = metric[s][m].omega;             
                e_gsm = exp(-0.5*gsm);
                e_rsm = exp(rsm);
                temp1 = trig.sin_theta[m];

                sum_rho += -e_gsm*trig.P_2n[m][0+1]*D2_metric[s][0+1].rho; 

                for(n=1; n<=LMAX; ++n) {

                    sum_rho += -e_gsm*trig.P_2n[m][n+1]*D2_metric[s][n+1].rho; 

                    if(m == MDIV) {             
                        sum_gama += -(2.0/PI)*e_gsm*D2_metric[s][n+1].gama;   
                        sum_omega += 0.5*e_rsm*e_gsm*D2_metric[s][n+1].omega; 
                    } else { 
                        sum_gama += -(2.0/PI)*e_gsm*(trig.sin_2n_1_theta[m][n]
                                    /((2.0*n-1.0)*temp1))*D2_metric[s][n+1].gama;   
                        sum_omega += -e_rsm*e_gsm*(trig.P1_2n_1[m][n+1]/(2.0*n
                                     *(2.0*n-1.0)*temp1))*D2_metric[s][n+1].omega;
                    }
                }
       
                metric[s][m].rho = rsm + control.cf*(sum_rho-rsm);
                metric[s][m].gama = gsm + control.cf*(sum_gama-gsm);
                metric[s][m].omega = omsm + control.cf*(sum_omega-omsm);

                sum_omega = 0.0;
                sum_rho = 0.0;
                sum_gama = 0.0; 
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &coeff_stop);
        printf("spin(), coeff: %ld\n", getElapsedTimeNs(coeff_start, coeff_stop));

        /* CHECK FOR DIVERGENCE */

        if(fabs(metric[2][1].rho) > 100.0 || fabs(metric[2][1].gama) > 300.0 || fabs(metric[2][1].omega) > 100.0) {
            a_check=200; 
            break;
        }


        /* TREAT SPHERICAL CASE */
      
        if(r_ratio == 1.0) {
            for(s=1; s<=SDIV; ++s) {
                for(m=1; m<=MDIV; ++m) {
                    metric[s][m].rho = metric[s][1].rho;
                    metric[s][m].gama = metric[s][1].gama;
                    metric[s][m].omega = 0.0;          
                }
            }
        }
      

        /* TREAT INFINITY WHEN SMAX=1.0 */

        if(SMAX == 1.0) {
            for(m=1; m<=MDIV; ++m) {
                metric[SDIV][m].rho = 0.0;
                metric[SDIV][m].gama = 0.0;
                metric[SDIV][m].omega = 0.0;
            }
        } 
      
        /* COMPUTE FIRST ORDER DERIVATIVES OF GAMA */ 
        struct timespec alpha_start, alpha_stop;
        clock_gettime(CLOCK_MONOTONIC, &alpha_start);

 
        for(s=1;s<=SDIV;++s) {
            for(m=1; m<=MDIV; ++m) {
                dgds[s][m] = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::gama, s, m);
                dgdm[s][m] = deriv_m<SDIV+1, MDIV+1>(metric, &Metric::gama, s, m);
            }
        }


        /* ALPHA */
 
        if(r_ratio == 1.0) {
            for(s=1; s<=SDIV; ++s) {
                for(m=1; m<=MDIV; ++m) {
                     da_dm[s][m] = 0.0; 
                }
            }
        } else {
            for(s=2; s<=SDIV; ++s) {
                for(m=1; m<=MDIV; ++m) {

                    da_dm[1][m] = 0.0; 
       
                    sgp = s_gp[s];
                    s1 = sgp*(1.0-sgp);
                    mum = mu[m]; 
                    m1 = 1.0-SQ(mum);
          
                    d_gama_s = dgds[s][m];
                    d_gama_m = dgdm[s][m];
                    d_rho_s = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::rho, s, m);
                    d_rho_m = deriv_m<SDIV+1, MDIV+1>(metric, &Metric::rho, s, m);
                    d_gama_sm = deriv_sm<SDIV+1, MDIV+1>(metric, &Metric::gama, s, m);
                    d_omega_s = deriv_s<SDIV+1, MDIV+1>(metric, &Metric::omega, s, m);
                    d_omega_m = deriv_m<SDIV+1, MDIV+1>(metric, &Metric::omega, s, m);
                    d_gama_ss = s1*deriv_s<double, SDIV+1, MDIV+1>(dgds, s, m)+(1.0-2.0*sgp)*d_gama_s;
                    d_gama_mm = m1*deriv_m<double, SDIV+1, MDIV+1>(dgdm, s, m)-2.0*mum*d_gama_m;  

                    temp1 = 2.0*SQ(sgp)*(sgp/(1.0-sgp))*m1*d_omega_s*d_omega_m
                            *(1.0+s1*d_gama_s) - (SQ(SQ(sgp)*d_omega_s) - 
                            SQ(sgp*d_omega_m/(1.0-sgp))*m1)*(-mum+m1*d_gama_m); 
  
                    temp2 = 1.0/(m1 *SQ(1.0+s1*d_gama_s) + SQ(-mum+m1*d_gama_m));

                    temp3 = s1*d_gama_ss + SQ(s1*d_gama_s);
  
                    temp4 = d_gama_m*(-mum+m1*d_gama_m);
   
                    temp5 = (SQ(s1*(d_rho_s+d_gama_s)) - m1*SQ(d_rho_m+d_gama_m))
                            *(-mum+m1*d_gama_m);

                    temp6 = s1*m1*(0.5*(d_rho_s+d_gama_s)* (d_rho_m+d_gama_m) 
                            + d_gama_sm + d_gama_s*d_gama_m)*(1.0 + s1*d_gama_s); 

                    temp7 = s1*mum*d_gama_s*(1.0+s1*d_gama_s);

                    temp8 = m1*exp(-2*metric[s][m].rho);
 
                    da_dm[s][m] = -0.5*(d_rho_m+d_gama_m) - temp2*(0.5*(temp3 - 
                                  d_gama_mm - temp4)*(-mum+m1*d_gama_m) + 0.25*temp5 
                                  - temp6 +temp7 + 0.25*temp8*temp1);     
                }
            }
        }

        for(s=1; s<=SDIV; ++s) {
            metric[s][1].alpha = 0.0;
            for(m=1; m<=MDIV-1; ++m) { 
                metric[s][m+1].alpha = metric[s][m].alpha+0.5*DM*(da_dm[s][m+1]+da_dm[s][m]);
            }
        } 
 
        for(s=1; s<=SDIV; ++s) {
            for(m=1; m<=MDIV; ++m) {     
                metric[s][m].alpha += -metric[s][MDIV].alpha+0.5*(metric[s][MDIV].gama-metric[s][MDIV].rho);

                if(metric[s][m].alpha >= 300.0) {
                    a_check = 200; 
                    break;
                }
                metric[s][m].omega /= r_e;
            }
        } 

        if(SMAX == 1.0) {
            for(m=1; m<=MDIV; ++m) {     
                metric[SDIV][m].alpha = 0.0;
            }
        }
        clock_gettime(CLOCK_MONOTONIC, &alpha_stop);
        printf("spin(), alpha: %ld\n", getElapsedTimeNs(alpha_start, alpha_stop));

        if(a_check == 200) {
            break;
        }

        dif = fabs(r_e_old-r_e)/r_e;
        n_of_it++;

    }   /* end while */



    /* COMPUTE OMEGA */  

    if(eos.isTabulatedEos()) { 
         Omega = Omega_h*C/(r_e*sqrt(KAPPA));
    } else {
         Omega = Omega_h/r_e;
    }

    /* UPDATE r_e_new */

    r_e_new = r_e;

}
