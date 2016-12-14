#include <string>

class EquationOfState {
    private:
        bool kIsTabulatedEos;
        double kGammaP;
        double _log_e_tab[201]; 
        double _log_p_tab[201];
        double _log_h_tab[201];
        double _log_n0_tab[201];
        int _num_tab; 
        void loadTabulatedEos(char eos_file[]);

    public:
        EquationOfState(); 
        EquationOfState(char eos_file[]); 
        EquationOfState(double gamma_p); 
        double e_of_rho0(double rho0) const;
        double e_at_p(double pp) const;         
        double e_at_p(double pp, int& n_nearest_pt) const;         
        double p_at_e(double ee, int& n_nearest_pt) const;        
        double p_at_h(double hh, int& n_nearest_pt) const;
        double h_at_p(double pp, int& n_nearest_pt) const;
        double n0_at_e(double ee, int& n_nearest_pt) const;
        bool isTabulatedEos() const;
        double getGammaP() const;
        int getNumTab() const;
};
