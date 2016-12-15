#include "RotatingNeutronStar.hh"


RotatingNeutonStar::RotatingNeutronStar(const char eos_file[], double central_energy_density) :
        eos(eos_file),
        e_center(central_energy_density) 
{
    initalizeRns();
}

RotatingNeutonStar::RotatingNeutronStar(double gamma_p, double central_energy_density) :
        eos(gamma_p),
        e_center(central_energy_density) 
{
    initializeRns();
}

void RotatingNeutronStar::initializeRns() {
    make_grid();
    make_center();
    make_surface();
    sphere();
}




