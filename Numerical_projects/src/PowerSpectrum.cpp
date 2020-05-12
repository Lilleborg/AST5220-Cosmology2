#include"PowerSpectrum.h"
#define W15 << std::setw(15) <<

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve()
{

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  // Vector k_array;
  // Vector log_k_array = log(k_array);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  // line_of_sight_integration(k_array);
  line_of_sight_integration_single();

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  // auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  // cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================
void PowerSpectrum::generate_bessel_function_splines()
{
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  
  // Argument for bessel functions
  const int n_arg = 2000;
  const double arg_min = 0;
  const double arg_max = 5000;
  Vector arg_array = Utils::linspace(arg_min,arg_max,n_arg);
  Vector j_ell_array(n_arg);
 
  // Loop over the different ells used
  for(int i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    // Loop over the range of argument values
    for (int i_arg = 0; i_arg < n_arg; i_arg++)
    {
      // Store the j_ell values to be used in the spline
      j_ell_array[i_arg] = Utils::j_ell(ell,arg_array[i_arg]);
    }
    // Make the j_ell_splines[i] spline
    std::string splinename("j_ell_spline_"+std::to_string(ell));
    j_ell_splines[i].create(arg_array,j_ell_array,splinename);
  }

  Utils::EndTiming("besselspline");
  // std::string data_path ("../data_testing/");
  // std::ofstream fp_bessel_spline (data_path + "test_bessel_spline.txt");
  // fp_bessel_spline W15 "arg";
  // for (int i = 0; i < ells.size(); i++)
  // {
  //   const int ell = ells[i];
  //   fp_bessel_spline W15 ell;
  // }
  // fp_bessel_spline << "\n";
  // for (int i_arg = 0; i_arg < n_arg; i_arg++)
  // {
  //   fp_bessel_spline W15 arg_array[i_arg];
  //   for (int i = 0; i < ells.size(); i++)
  //   {
  //     fp_bessel_spline W15 j_ell_splines[i](arg_array[i_arg]);
  //   }
  //   fp_bessel_spline << "\n";
  // }
  // fp_bessel_spline.close();
}

// int PowerSpectrum::rhs_theta_ell_ode(double x, double k, )
//====================================================
// Do the line of sight integration for a single
// source function
//====================================================
void PowerSpectrum::line_of_sight_integration_single()
{
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  // Set up initial condition for Theta_ell equal zero for x_start
  Vector Theta_ell_IC{0};
  // Set up ODESolver object to be used for the integration
  double hstart = 1e-4, abserr = 1e-10, relerr = 1e-10;
  ODESolver ODE_theta_ell(hstart,abserr,relerr);

  // Loop over the different ks
  for(int ik = 0; ik < k_array.size(); ik++)
  {
    const double k = k_array[ik];
    // Loop over the different ells
    for (int iell = 0; iell < ells.size(); iell++)
    {
      // Set up ODE rhs for this k and ell
      ODEFunction dtheta_ell_dx = [&](double x, const double *y, double *dydx)
      {
        const double arg = k*(cosmo->eta_of_x_spline(0.0)-cosmo->eta_of_x_spline(x));
        dydx[0] = pert->get_Source_T(x,k)*j_ell_splines[iell](arg);
        return GSL_SUCCESS;
      };
      // Solve the ODE
      ODE_theta_ell.solve(dtheta_ell_dx,x_array,Theta_ell_IC);
      // Extract result today
      result[iell][ik] = ODE_theta_ell.get_final_data_by_component(0);
    }
  }
  Utils::EndTiming("lineofsight");
  
  thetaT_ell_of_k_spline2D.create(ells,k_array,result,"thetaT_ell_of_k_spline2D");
}

//====================================================
// Do the line of sight integration
//====================================================
// void PowerSpectrum::line_of_sight_integration(Vector & k_array)
// {
//   const int n_k        = k_array.size();
//   const int n          = 100;
//   const int nells      = ells.size();
  
//   // Make storage for the splines we are to create
//   thetaT_ell_of_k_spline = std::vector<Spline>(nells);

//   //============================================================================
//   // TODO: Solve for Theta_ell(k) and spline the result
//   //============================================================================

//   // Make a function returning the source function
//   std::function<double(double,double)> source_function_T = [&](double x, double k){
//     return pert->get_Source_T(x,k);
//   };

//   // Do the line of sight integration
//   Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

// }

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  // ...
  // ...
  // ...
  // ...

  Vector result;

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

