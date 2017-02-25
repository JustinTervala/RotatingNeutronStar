#include "GridTrig.hh"
#include "EquationOfState.hh"
#include "equil_util.hh"


//TODO: Use static polymorphism to do a polymorphic and a tabulated RNS
class RotatingNeutronStar {
    private:
        EquationOfState eos;
        double e_center, p_center, h_center;
        double e_surface,  p_surface;
        double enthalpy_min;
        double r_e;
        double Omega;
        double Omega_K;
        double Mass; 
        double Mass_0;
        double J;
        double r_ratio;
        double ang_mom;
        double R_e;
        double accuracy;
        double cf;
        double s_gp[SDIV+1];
        double mu[MDIV+1];
        GridTrig trig;
        matrix<Metric, SDIV+1, MDIV+1> metric;
        matrix<double, SDIV+1, MDIV+1> energy;
        matrix<double, SDIV+1, MDIV+1> pressure;
        matrix<double, SDIV+1, MDIV+1> enthalpy;
        matrix<double, SDIV+1, MDIV+1> velocity_sq;
        std::array<double, SDIV+1> v_plus;
        std::array<double, SDIV+1> v_minus;
        void initializeRns();
        void make_grid();
        void make_center(); 

        void TOV(int    i_check, 
                 double r_is_gp[RDIV+1], 
                 double lambda_gp[RDIV+1], 
                 double nu_gp[RDIV+1], 
                 double &r_is_final, 
                 double &r_final, 
                 double &m_final);
        
        double dm_dr_is(double r_is, 
                        double r, 
                        double m, 
                        double p, 
                        int    &n_nearest_pt);

        double dp_dr_is(double r_is, 
                        double r, 
                        double m, 
                        double p,
                        int    &n_nearest_pt);

        double dr_dr_is(double r_is, double r, double m);
    public:
        RotatingNeutronStar(char eos_file[], double central_energy_density);
        RotatingNeutronStar(double gamma_p, double central_energy_density);
        void sphere();
        void spin();
        void mass_radius();
        void recompute(double r_ratio_);
        void print_state() const;
        void setAccuracy(double accuracy_);
        void setCf(double cf_); 
        void setRRatio(double r_ratio_);
        double getOmega() const;
        double getOmegaK() const;
};
