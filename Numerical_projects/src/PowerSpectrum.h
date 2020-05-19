#ifndef _POWERSPECTRUM_HEADER
#define _POWERSPECTRUM_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <functional>
#include <utility> 
#include <fstream> 
#include <algorithm>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class PowerSpectrum {
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
    Perturbations *pert        = nullptr;

    // Parameters defining the primordial power-spectrum
    double A_s        = 2e-9;
    double n_s        = 0.96;
    double kpivot_mpc = 0.05;

    // Bool to determine if all the different terms in the Source function is computed separately
    bool treat_source_components = false;

    // The k-values we compute Theta_ell(k) etc. for
    const int n_k      = 2000;
    const double k_min = Constants.k_min;
    const double k_max = Constants.k_max;
    Vector k_array     = Utils::logspace(log(k_min),log(k_max),n_k);
    Vector log_k_array = log(k_array);

    // The x vector
    int n_x;
    const double x_min = Constants.x_start;
    const double x_max = 0;

    // The ells's we will compute Theta_ell and Cell for
    Vector ells{ 
        2,    3,    4,    5,    6,    7,    8,    10,   12,   15,   
        20,   25,   30,   40,   50,   60,   70,   80,   90,   100,  
        120,  140,  160,  180,  200,  225,  250,  275,  300,  350,  
        400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  
        900,  950,  1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 
        1900, 1950, 2000};
   
    //=====================================================================
    // [1] Create bessel function splines needed for the LOS integration
    //=====================================================================

    // Splines of bessel-functions for each value of ell in the array above
    std::vector<Spline> j_ell_splines;
    
    // Generate splines of bessel-functions for each ell needed
    // to do the LOS integration
    void generate_bessel_function_splines();
    
    //=====================================================================
    // [2] Do the line of sight integration and spline the result
    //=====================================================================
    
    // Do LOS integration for all ells and all k's in the given k_array
    // and for all the source functions (temperature, polarization, ...)
    void line_of_sight_integration(bool solve_source_components, const Vector x);
  
    // Do the line of sight integration for a single quantity
    // for all ells by providing a source_function(x,k) (can be temp, pol, ...)
    Vector2D line_of_sight_integration_single(std::function<double(double,double)> &source_function, const Vector x);
    
    // Splines of the result of the LOS integration
    // Theta_ell(k) and ThetaE_ell(k) for polarization
    Spline2D thetaT_ell_of_k_spline2D{"thetaT_ell_of_k_spline2D"};
    Spline2D thetaSW_ell_of_k_spline2D{"thetaSW_ell_of_k_spline2D"};
    Spline2D thetaISW_ell_of_k_spline2D{"thetaISW_ell_of_k_spline2D"};
    Spline2D thetaDoppler_ell_of_k_spline2D{"thetaDoppler_ell_of_k_spline2D"};
    Spline2D thetaQuad_ell_of_k_spline2D{"thetaQuad_ell_of_k_spline2D"};
    Spline2D thetag_tilde_ell_of_k_spline2D{"thetag_tilde_ell_of_k_spline2D"};
    
    //=====================================================================
    // [3] Integrate to get power-spectrum
    //=====================================================================
    
    // General method to solve for Cells (allowing for cross-correlations)
    // For auto spectrum (C_TT) then call with f_ell = g_ell = theta_ell
    // For polarization C_TE call with f_ell = theta_ell and g_ell = thetaE_ell
    Vector solve_for_cell(
        Vector & logk_array,
        Spline2D & f_ell, 
        Spline2D & g_ell);

    // Splines with the power-spectra
    Spline cell_TT_spline{"cell_TT_spline"};
    Spline cell_SW_spline{"cell_SW_spline"};
    Spline cell_ISW_spline{"cell_ISW_spline"};
    Spline cell_Doppler_spline{"cell_Doppler_spline"};
    Spline cell_Quad_spline{"cell_Quad_spline"};
    Spline cell_g_tilde_spline{"cell_g_tilde_spline"};
    Spline cell_TE_spline{"cell_TE_spline"};
    Spline cell_EE_spline{"cell_EE_spline"};

  public:

    // Constructors
    PowerSpectrum() = delete;
    PowerSpectrum(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec, 
        Perturbations *pert);
    
    // Do all the solving: bessel functions, LOS integration and then compute Cells
    void solve(bool solve_source_components = false);

    // The dimensionless primordial power-spectrum Delta = 2pi^2/k^3 P(k)
    double primordial_power_spectrum(const double k) const;

    // Get P(k,x) for a given x in units of (Mpc)^3
    double get_component_power_spectrum(const std::string component, const double x, const double k_mpc) const;

    // Gauge invariant density perturbation for "component"
    double get_invariant_density(const std::string component, const double x, const double k) const;

    // Get the quantities we have computed
    double get_cell_TT(const double ell) const;
    double get_cell_TE(const double ell) const;
    double get_cell_EE(const double ell) const;

    // Output Cells in units of l(l+1)/2pi (muK)^2
    void output(std::string filename) const;
    void output_component_power_spectrum(std::vector<std::string> components, std::string filename) const;
};

#endif
