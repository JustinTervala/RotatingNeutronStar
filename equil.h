#include <time.h>
#include <array>
#include "equil_util.h"
#include "EquationOfState.hh"

long getElapsedTimeNs(struct timespec start, struct timespec stop);

void make_grid(double s_gp[SDIV+1], 
               double mu[MDIV+1]);                        

void make_center(const EquationOfState& eos, 
                 double e_center,
                 double &p_center, 
                 double &h_center);

void mass_radius(double s_gp[SDIV+1],
                 double mu[MDIV+1],
                 const EquationOfState& eos, 
                 std::array<std::array<Metric, MDIV+1>, SDIV+1>& metric,
                 std::array<std::array<double, MDIV+1>, SDIV+1>& energy,
                 std::array<std::array<double, MDIV+1>, SDIV+1>& pressure,
                 std::array<std::array<double, MDIV+1>, SDIV+1>& enthalpy,
                 std::array<std::array<double, MDIV+1>, SDIV+1>& velocity_sq,
                 double r_ratio,
                 double e_surface,
                 double r_e,
                 double Omega,
                 double &Mass, 
                 double &Mass_0,
                 double &ang_mom,
                 double &R_e,
                 std::array<double, SDIV+1>& v_plus,
                 std::array<double, SDIV+1>& v_minus,
                 double &Omega_K);

double dm_dr_is(double r_is, 
                double r, 
                double m, 
                double p, 
                double e_center, 
                double p_surface,
                const EquationOfState& eos, 
                int    &n_nearest_pt);

double dp_dr_is(double r_is, 
                double r, 
                double m, 
                double p,
                double e_center, 
                double p_surface,
                const EquationOfState& eos,
                int    &n_nearest_pt);

double dr_dr_is(double r_is, double r, double m);

void TOV(int    i_check, 
         double e_center,
         double p_center,
         double p_surface,
         double e_surface,
         const EquationOfState& eos, 
         double r_is_gp[RDIV+1], 
         double lambda_gp[RDIV+1], 
         double nu_gp[RDIV+1], 
         double &r_is_final, 
         double &r_final, 
         double &m_final);

void sphere(double s_gp[SDIV+1],
            const EquationOfState& eos, 
            double e_center,
            double p_center, 
            double h_center,
            double p_surface,
            double e_surface,
            std::array<std::array<Metric, MDIV+1>, SDIV+1>& metric,
            double &r_e);


void spin(double s_gp[SDIV+1],
          double mu[MDIV+1],
          const EquationOfState& eos, 
          double h_center,
          double enthalpy_min,
          std::array<std::array<Metric, MDIV+1>, SDIV+1>& metric,
          std::array<std::array<double, MDIV+1>, SDIV+1>& energy,
          std::array<std::array<double, MDIV+1>, SDIV+1>& pressure,
          std::array<std::array<double, MDIV+1>, SDIV+1>& enthalpy,
          std::array<std::array<double, MDIV+1>, SDIV+1>& velocity_sq,
          int    a_check, 
          double accuracy,
          double cf,
          double r_ratio,
          double &r_e_new,
          double &Omega) ;
 
void compute_f2n(const double s_gp[SDIV+1],
                 double (&f2n)[LMAX+2][SDIV+1]);

void compute_f_rho_gamma(const double s_gp[SDIV+1],
                         const double (&f2n)[LMAX+2][SDIV+1],
                         float (&f_rho)[SDIV+1][LMAX+2][SDIV+1],
                         float (&f_gama)[SDIV+1][LMAX+2][SDIV+1]);

void compute_trig(const double mu[MDIV+1],
                  double theta[MDIV+1],
                  double sin_theta[MDIV+1],
                  double (&P_2n)[MDIV+1][LMAX+2],
                  double (&P_2n_t)[LMAX+2][MDIV+1],
                  double (&P1_2n_1)[MDIV+1][LMAX+2],
                  double (&P1_2n_1_t)[LMAX+2][MDIV+1],
                  double (&sin_2n_1_theta)[MDIV+1][LMAX+1],
                  double (&sin_2n_1_theta_t)[LMAX+1][MDIV+1]);
                   
