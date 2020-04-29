#ifndef _PERTURBATIONS_HEADER
#define _PERTURBATIONS_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <vector>
#include <fstream>
#include <algorithm>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class Perturbations{
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
   
    // The scales we integrate over
    const int n_k        = 150; //use only one k-value to begin with   was: 100;
    const double k_min   = Constants.k_min;
    const double k_max   = Constants.k_max;
    // Set up logarithmic (using base e as for logarithmic scale factor) 
    // spaced k-values using own logspace from Utils
    Vector k_array = Utils::logspace(log(k_min),log(k_max),n_k);
    
    // Start and end of the time-integration
    const int n_x        = 2000;
    const double x_start = Constants.x_start;
    const double x_end   = Constants.x_end;
    // Set up logarithmic scale factor
    Vector x_array_full = Utils::linspace(x_start,x_end,n_x);
    
    // Below is a full list of splines you probably need, 
    // but you only need to make the splines you will need

    // Splines of scalar perturbations quantities
    Spline2D delta_cdm_spline{"delta_cdm_spline"};
    Spline2D delta_b_spline{"delta_b_spline"};
    Spline2D v_cdm_spline{"v_cdm_spline"};
    Spline2D v_b_spline{"v_b_spline"};
    Spline2D dv_b_dx_spline{"dv_b_dx_spline"};
    Spline2D Pi_spline{"Pi_spline"};
    Spline2D Psi_spline{"Psi_spline"};
    Spline2D Phi_spline{"Phi_spline"};
    Spline2D dPhi_dx_spline{"dPhi_dx_spline"};
   
    // Splines of source functions (ST for temperature; SE for polarization)
    Spline2D ST_spline{"ST_spline"};
    // Spline2D SE_spline{"SE_spline"};
    
    // Splines of multipole quantities
    Spline2D Theta0_spline{"Theta0_spline"};
    Spline2D Theta1_spline{"Theta1_spline"};
    Spline2D Theta2_spline{"Theta2_spline"};
    std::vector<Spline2D> vector_of_Theta_splines;
    Spline2D dTheta0_dx_spline{"Theta0_spline"};
    Spline2D dTheta1_dx_spline{"Theta1_spline"};
    std::vector<Spline2D> vector_of_dTheta_dx_splines;

    
    //==========================================================
    // [1] Tight coupling ODE system
    //==========================================================

    // Set the initial conditions at the start (which is in tight coupling)
    Vector set_ic(const double x, const double k) const;
    
    // Right hand side of the ODE in the tight coupling regime  (A vector from slides, only few multipoles of theta)
    int rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx);
    
    // Compute the time when tight coupling ends  (what x does tight coupling end for each k)
    std::pair<double,int> get_tight_coupling_time_and_index(const double k) const;
    
    //==========================================================
    // [2] The full ODE system 
    //==========================================================
    
    // Set initial condition after tight coupling
    Vector set_ic_after_tight_coupling(const Vector &y_tight_coupling, const double x, const double k) const;

    // Right hand side of the ODE in the full regime (A vector from slides, all multiploles of theta)
    int rhs_full_ode(double x, double k, const double *y, double *dydx);
    
    //==========================================================
    // [3] Integrate the full system
    //==========================================================
    
    // Integrate perturbations and spline the result
    void integrate_perturbations();
    
    //==========================================================
    // [4] Compute source functions from the result
    //==========================================================
    
    // Compute source functions and spline the result
    void compute_source_functions();

  public:

    // Constructors
    Perturbations() = default;
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec); 

    // Set up x_array with higher resolution during recombination, following Callin
    Vector set_up_x_array_resolution() const;

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output info to file
    void output(const double k, const std::string filename) const;

    // Get the quantities we have integrated
    double get_delta_cdm(const double x, const double k) const;
    double get_delta_b(const double x, const double k) const;
    double get_v_cdm(const double x, const double k) const;
    double get_v_b(const double x, const double k) const;
    double get_Phi(const double x, const double k) const;
    double get_Psi(const double x, const double k) const;
    double get_Pi(const double x, const double k) const;
    double get_Theta(const double x, const double k, const int ell) const;
    // double get_Theta_p(const double x, const double k, const int ell) const;
    // double get_Nu(const double x, const double k, const int ell) const;
    double get_Source_T(const double x, const double k) const;
    // double get_Source_E(const double x, const double k) const;
};

#endif
