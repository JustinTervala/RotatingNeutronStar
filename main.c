/************************************************************************** 
*                            MAIN.C                                       * 
*                                                                         *
*       This is a sample program which uses the routines 
*       written by Nikolaos Stergioulas. 
*
* This sample program is meant merely to serve as an 
* example of how to use the functions! 
*
* In this example, the user specifies the star's equation of state,
* and the central energy density. The program computes models with
* increasing angular velocity until the star is spinning with 
* the same angular velocity as a particle orbiting the star
* at its equator. The last star computed will be spinning with
* the maximum allowed angular velocity for a star with the given
* central energy density. 
*
* Note that there is no guarantee that any of the intermediate
* stars are stable to radial perturbations. Also there is no
* guarantee that given any energy density, there will be a 
* stable rotating solution. As a rule of thumb, the highest 
* energy density you should use is the value which gives the
* maximum mass nonrotating star. 
*
* It would be a good idea to read some of the papers on rapidly
* rotating neutron stars (such as given on the rns homepage)
* before embarking on a study of rotating neutron stars
* using these routines. 
* 
* For example, a reasonable model would use the file eosA,
* central energy density = 10^{15} g/cm^3 
* To specify these parameters, run the executable:
 
kepler -f eosA -e 1e15 

* This program (compiled on the "standard" setting) 
* requires about 2.7 MBytes of RAM and took about 2 minutes to run.
*                                                                         *
**************************************************************************/
#include "consts.h"
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <array>
#include "GridTrig.hh"
#include "RotatingNeutronStar.hh"

/*************************************************************************/
/* Main program.                                                         */
/*************************************************************************/
int main(int argc,                    /* Number of command line arguments */ 
         char **argv) {               /* Command line arguments */

 /* EQUILIBRIUM VARIABLES */

 int n_tab;                     /* Number of points in EOS file */
       
 struct timespec spin_start, spin_stop, mr_start, mr_stop;

 int i, j;

 double log_e_tab[201],               /* energy density/c^2 in tabulated EOS */
        log_p_tab[201],               /* pressure in tabulated EOS */
        log_h_tab[201],               /* enthalpy in EOS file */
        log_n0_tab[201],              /* number density in EOS file */  
        e_center,                     /* central en. density */
        p_center,                     /* central pressure */
        h_center,                     /* central enthalpy */
        n_P,                          /* Polytropic index N */
        Gamma_P,                      /* Gamma for polytropic EOS */     
        r_ratio,                      /* axis ratio */
        s_gp[SDIV+1],                 /* s grid points */
        mu[MDIV+1],                   /* \mu grid points */
        Mass,                         /* Gravitational mass */
        e_surface,                    /* surface en. density */ 
        p_surface,                    /* surface pressure */
        enthalpy_min,                 /* minimum enthalpy in EOS */
        Mass_0,                       /* Baryon Mass */
        Omega,                        /* Angular Velocity */
        J,                            /* Angular Momentum */
        R_e,                          /* Circumferential radius at equator */
       //*v_plus,                       /* vel. of co-rot. particle wrt ZAMO */
       //*v_minus,                      /* vel. of counter-rot. ... */
        Omega_K,                      /* Keplerian velocity of particle orbiting at equator */
        r_e;                          /* coord. radius at equator     */


    double cf=1, /* convergence factor */
           accuracy;

    /* Definitions used to find the interval containing the 
       correct spin frequency */

    double dr,
           diff_Omega,
           old_diff_Omega,
           a_check;
   
    /* Definitions used in the Ridder zero finding method */
 
    float ans, fh, fl, fm, fnew, sroot, xh, xl, xm, xnew, xacc;

    char eos_file[80] = "no EOS file specified";   /* EOS file name */
    char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
 

    /* READ IN THE COMMAND LINE OPTIONS */

    for(i=1; i<argc; i++) { 
        if(argv[i][0] == '-'){
            switch(argv[i][1]){

                case 'q':
                    /* CHOOSE THE EOS TYPE: EITHER "tab" or "poly"(default is tab) */
                    sscanf(argv[i+1], "%s", eos_type);
                    break;
               
                case 'N': 
                    /* IF A POLYTROPIC EOS WAS CHOSEN, CHOOSE THE POLYTROPIC INDEX "N" */
                    sscanf(argv[i+1], "%lf", &n_P);
                    Gamma_P = 1.0+1.0/n_P;
                    break;               

                case 'f':
                    /* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE NAME OF THE FILE */
                    sscanf(argv[i+1], "%s", eos_file);
                    break;

                case 'e':
                    /* CHOOSE THE CENTRAL ENERGY DENSITY OF THE NEUTRON STAR (IN g/cm^3) */
                    sscanf(argv[i+1], "%lf", &e_center);
                    if(strcmp(eos_type, "tab") == 0) {
                        e_center *= C*C*KSCALE;
                    }
                    break;
                    
                case 'h': 
                    fprintf(stderr, "\n");
                    fprintf(stderr, "Quick help:\n");
                    fprintf(stderr, "\n");
                    fprintf(stderr, "  -q EOS type (tab)\n"); 
                    fprintf(stderr, "     tab : tabulated \n");
                    fprintf(stderr, "     poly : analytic polytropic \n");           
                    fprintf(stderr, "  -N polytropic index (P=K*e^(1+1/N))\n");  
                    fprintf(stderr, "  -f EOS file \n");
                    fprintf(stderr, "  -e central energy density in gr/cm^3\n");
                    fprintf(stderr, "  -h this menu\n");
                    fprintf(stderr, "\n");
                    exit(1);
                    break;
     
            }
        }
    }

    /* PRINT THE HEADER */

    if(strcmp(eos_type, "tab") == 0) {
        printf("%s,  MDIVxSDIV=%dx%d\n", eos_file, MDIV, SDIV);
    } else {
        printf("polytrope with N=%f, MDIVxSDIV=%dx%d\n", n_P, MDIV, SDIV);
    }
    printf("ratio\te_15\tM\tM_0\tr_star\tspin\tOmega_K\tI\tJ/M^2\n");
    if(strcmp(eos_type, "tab") == 0) {
        printf("\tg/cm^3\tsun\tsun\tkm\ts-1\ts-1\tg cm^2\t\n");
    }
    printf("\n");

    RotatingNeutronStar rns = (strcmp(eos_type, "tab") == 0) ? RotatingNeutronStar(eos_file, e_center) : RotatingNeutronStar(Gamma_P, e_center);
    /* ALLLOCATE MEMORY */

    /* set program defaults */
    rns.setAccuracy(1e-5);  
    rns.setCf(1.0); 
    xacc = 1e-4;  
 
    r_ratio = 1.0; 
    rns.recompute(r_ratio);    
    
    /* PRINT OUT INFORMATION ABOUT THE STELLAR MODEL */
    rns.print_state();
    dr = 0.05;
 

    /* THIS LOOP STARTS WITH A NON-ROTATING STAR AND INCREASES
       THE STAR'S OBLATENESS (BY DECREASING R_RATIO) AND 
       THEN CALCULATES THE STAR'S ANGULAR VELOCITY. ONCE THE
       COMPUTED VALUE OF ANGULAR VELOCITY IS LARGER THAN 
       THE ANGULAR VELOCITY OF A PARTICLE ORBITING THE STAR
       AT THE EQUATOR, (Omega_K), THE LOOP STOPS */
    
    diff_Omega = rns.getOmegaK() - rns.getOmega();
    old_diff_Omega = diff_Omega;
  
    while(diff_Omega > 0) {
        /* Find the interval of r_ratio where the star has the
           correct angular velocity    */
        r_ratio -= dr;
        rns.recompute(r_ratio);    

        rns.print_state();
   
        old_diff_Omega = diff_Omega;
        diff_Omega = rns.getOmegaK() - rns.getOmega();
    } 

    /* The correct star lies between r_ratio and r_ratio + dr */
    xl = r_ratio;
    xh = r_ratio + dr;
    fl = diff_Omega;
    fh = old_diff_Omega;

    /* Use Ridder's method to find the correct star (Taken from 
       Numerical Recipes)    */

    if((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
        ans = -1.11e30;
        for (j=1; j<=60; j++) {
            xm = 0.5*(xl+xh);
            r_ratio = xm;
            rns.recompute(r_ratio);    
            
            fm = rns.getOmegaK() - rns.getOmega();
            sroot = sqrt(fm*fm-fl*fh);
            
            if(sroot == 0.0) {
                r_ratio = ans;
                break;
            }
           
            xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/sroot);
            if(fabs(xnew-ans) <= xacc) {
                r_ratio = ans;
                break;
            }
            ans = xnew;
            r_ratio = ans;
            rns.recompute(r_ratio);    
            rns.print_state();

            fnew =  rns.getOmegaK() - rns.getOmega();
            if(fnew == 0.0){
                r_ratio = ans;
                rns.setRRatio(r_ratio);
                break;
            }
           
            if(sgn<float>(fm, fnew) != fm) {
                xl = xm;
                fl = fm;
                xh = ans;
                fh = fnew;
            } else if(sgn<float>(fl, fnew) != fl) {
                xh = ans;
                fh = fnew;
            } else if(sgn<float>(fh, fnew) != fh) {
                xl = ans;
                fl = fnew;
            } else {
                print_error("never get here.");
            }
            if(fabs(xh-xl) <= xacc) {
                r_ratio = ans;
                rns.setRRatio(r_ratio);
                break;
            }
        }
    } else {
        if(fh == 0.0){
            r_ratio +=dr;
           rns.setRRatio(r_ratio);
        }
        print_error("root must be bracketed in zriddr.");
    }
 
    /* THE RIDDER ZERO-FINDING ROUTINE HAS FOUND THE VALUE
       OF R_RATIO WHICH GIVES THE DESIRED STAR. */

    return 0;

}









