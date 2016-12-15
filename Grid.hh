#pragma once

#include "consts.h"

struct Grid {
    double s_gp[SDIV+1]; 
    double mu[MDIV+1];
    Grid() {
        for(int s = 1; s <= SDIV; ++s) {
            s_gp[s] = SMAX*(s-1.0)/(SDIV-1.0);
        }
        /* s_gp[1] = 0.0     corresponds to the center of the star
           s_gp[SDIV] = SMAX corresponds to infinity */
        for(int m = 1; m<= MDIV; ++m) { 
             mu[m] = (m-1.0)/(MDIV-1.0);
        }
        /* mu[1] = 0.0    corresponds to the plane of the equator 
           mu[MDIV] = 1.0 corresponds to the axis of symmetry */

        /* s_gp[0] and mu[0] are not used by the program */
    }
    ~Grid() {}
};


