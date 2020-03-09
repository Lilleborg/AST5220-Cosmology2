#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "BackgroundCosmology.h"

using Vector = std::vector<double>;
using Doublepair = std::pair<double,double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;
 
    // The start and end points for recombination arrays (can be modified)
    const double x_start  = Constants.x_start;
    const double x_end    = Constants.x_end;
    
    // Numbers of points used to solved for recombination (Xe and ne)
    const int npts_rec_arrays = 1e+5;
    // Numbers of points used to solved for optical depth (tau)
    const int npts_tau        = 1e+3;
  
    // Xe for when to switch between Saha and Peebles
    const double Xe_saha_limit = 0.99;

    // Compute Xe from the Saha equation
    Doublepair electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe using Saha and Peebles
    void solve_number_density_electrons();
    
    // Solve for optical depth using the ODESolver, compute visibility function
    void solve_for_optical_depth_tau();

    // Splines contained in this class
    Spline log_Xe_of_x_spline{"log_Xe"};
    Spline log_ne_of_x_spline{"log_ne"};
    Spline tau_of_x_spline{"tau"};
    Spline tau_deriv_of_x_spline{"tau_deriv"};
    Spline g_tilde_of_x_spline{"g"};  
    Spline g_tilde_deriv_of_x_spline{"g"};  

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp);

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Print results of last scattering from recombination class
    void print_time_results() const;

    // Output some data to file
    void output(const std::string filename) const;

    // Get functions that we must implement
    double tau_of_x(double x) const;
    double dtaudx_of_x(double x) const;
    double ddtauddx_of_x(double x) const;
    double g_tilde_of_x(double x) const;
    double dgdx_tilde_of_x(double x) const;
    double ddgddx_tilde_of_x(double x) const;
    double Xe_of_x(double x) const;
    double ne_of_x(double x) const;
    double get_Yp() const;

    double get_number_density_baryons(double x) const;
    Vector get_time_results() const;
};

#endif
