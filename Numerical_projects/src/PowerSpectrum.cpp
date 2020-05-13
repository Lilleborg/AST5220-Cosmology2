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
  // Make splines for j_ell. 
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // Line of sight integration to get Theta_ell(k)
  //=========================================================================
  line_of_sight_integration();

  //=========================================================================
  // Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  //=========================================================================
  
  Vector cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline2D, thetaT_ell_of_k_spline2D);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  Vector cell_SW = solve_for_cell(log_k_array, thetaSW_ell_of_k_spline2D, thetaSW_ell_of_k_spline2D);
  cell_SW_spline.create(ells, cell_SW, "Cell_SW_of_ell");
  
  Vector cell_ISW = solve_for_cell(log_k_array, thetaISW_ell_of_k_spline2D, thetaISW_ell_of_k_spline2D);
  cell_ISW_spline.create(ells, cell_ISW, "Cell_ISW_of_ell");
  
  Vector cell_Doppler = solve_for_cell(log_k_array, thetaDoppler_ell_of_k_spline2D, thetaDoppler_ell_of_k_spline2D);
  cell_Doppler_spline.create(ells, cell_Doppler, "Cell_Doppler_of_ell");
  
  Vector cell_Quad = solve_for_cell(log_k_array, thetaQuad_ell_of_k_spline2D, thetaQuad_ell_of_k_spline2D);
  cell_Quad_spline.create(ells, cell_Quad, "Cell_Quad_of_ell");
  
  Vector cell_g_tilde = solve_for_cell(log_k_array, thetag_tilde_ell_of_k_spline2D, thetag_tilde_ell_of_k_spline2D);
  cell_g_tilde_spline.create(ells, cell_g_tilde, "Cell_g_tilde_of_ell");
  
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
  const int n_arg = 10000;
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

//====================================================
  // Do the line of sight integration for a single
  // source function
  //====================================================
Vector2D PowerSpectrum::line_of_sight_integration_single(std::function<double(double,double)> &source_function)
{
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  // Set up initial condition for Theta_ell equal zero for x_start
  Vector Theta_ell_IC{0};
  // Set up ODESolver object to be used for the integration
  double hstart = 1e-3, abserr = 1e-10, relerr = 1e-10;
  ODESolver ODE_theta_ell(hstart,abserr,relerr);

  // Loop over the different ks
  for(int ik = 0; ik < k_array.size(); ik++)
  {
    const double k = k_array[ik];
    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k )
    {
      printf("Progress pert integration: %3d%%, k-value per Mpc: %.3e\n",(100*ik+100)/n_k,k/Constants.Mpc);
      if(ik == n_k-1) std::cout << std::endl;
    }
    // Loop over the different ells
    for (int iell = 0; iell < ells.size(); iell++)
    {
      // Set up ODE rhs for this k and ell
      ODEFunction dtheta_ell_dx = [&](double x, const double *y, double *dydx)
      {
        const double arg = k*(cosmo->eta_of_x_spline(0.0)-cosmo->eta_of_x_spline(x));
        dydx[0] = source_function(x,k)*j_ell_splines[iell](arg);
        return GSL_SUCCESS;
      };
      // Solve the ODE
      ODE_theta_ell.solve(dtheta_ell_dx,x_array,Theta_ell_IC);
      // Extract result today
      result[iell][ik] = ODE_theta_ell.get_final_data_by_component(0);
    }
  }
  Utils::EndTiming("lineofsight");
  return result;
}

// ====================================================
// Do the line of sight integration
// ====================================================
void PowerSpectrum::line_of_sight_integration()
{
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();

  std::cout << "Solving with the full temperature source function\n";
  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };
  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(source_function_T);
  // Create thetaT 2D spline
  thetaT_ell_of_k_spline2D.create(ells,k_array,thetaT_ell_of_k,"thetaT_ell_of_k_spline2D");

  std::cout << "Solving with the SW source function\n";
  // Make a function returning the source function
  std::function<double(double,double)> source_function_SW = [&](double x, double k){
    return pert->get_Source_SW(x,k);
  };
  // Do the line of sight integration
  Vector2D thetaSW_ell_of_k = line_of_sight_integration_single(source_function_SW);
  // Create thetaSW 2D spline
  thetaSW_ell_of_k_spline2D.create(ells,k_array,thetaSW_ell_of_k,"thetaSW_ell_of_k_spline2D");

  std::cout << "Solving with the ISW source function\n";
  // Make a function returning the source function
  std::function<double(double,double)> source_function_ISW = [&](double x, double k){
    return pert->get_Source_ISW(x,k);
  };
  // Do the line of sight integration
  Vector2D thetaISW_ell_of_k = line_of_sight_integration_single(source_function_ISW);
  // Create thetaISW 2D spline
  thetaISW_ell_of_k_spline2D.create(ells,k_array,thetaISW_ell_of_k,"thetaISW_ell_of_k_spline2D");

  std::cout << "Solving with the Doppler source function\n";
  // Make a function returning the source function
  std::function<double(double,double)> source_function_Doppler = [&](double x, double k){
    return pert->get_Source_Doppler(x,k);
  };
  // Do the line of sight integration
  Vector2D thetaDoppler_ell_of_k = line_of_sight_integration_single(source_function_Doppler);
  // Create thetaDoppler 2D spline
  thetaDoppler_ell_of_k_spline2D.create(ells,k_array,thetaDoppler_ell_of_k,"thetaDoppler_ell_of_k_spline2D");

  std::cout << "Solving with the Quad source function\n";
  // Make a function returning the source function
  std::function<double(double,double)> source_function_Quad = [&](double x, double k){
    return pert->get_Source_Quad(x,k);
  };
  // Do the line of sight integration
  Vector2D thetaQuad_ell_of_k = line_of_sight_integration_single(source_function_Quad);
  // Create thetaQuad 2D spline
  thetaQuad_ell_of_k_spline2D.create(ells,k_array,thetaQuad_ell_of_k,"thetaQuad_ell_of_k_spline2D");

  std::cout << "Solving with the g_tilde source function\n";
  // Make a function returning the source function
  std::function<double(double,double)> source_function_g_tilde = [&](double x, double k){
    return rec->g_tilde_of_x(x);
  };
  // Do the line of sight integration
  Vector2D thetag_tilde_ell_of_k = line_of_sight_integration_single(source_function_g_tilde);
  // Create thetag_tilde 2D spline
  thetag_tilde_ell_of_k_spline2D.create(ells,k_array,thetag_tilde_ell_of_k,"thetag_tilde_ell_of_k_spline2D");
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    Spline2D & f_ell_spline,
    Spline2D & g_ell_spline)
{
  Utils::StartTiming("C_ell integration");
  const int nells = ells.size();
  // Storeage for result
  Vector result(nells);
  // Set up initial condition for C_ell equal zero for k_min
  Vector C_ell_IC{0};
  // Set up ODESolver object to be used for the integration
  double hstart = 1e-3, abserr = 1e-10, relerr = 1e-10;
  ODESolver ODE_C_ell(hstart,abserr,relerr);

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  for (int iell = 0; iell < nells; iell++)
  {
    const double ell = ells[iell];
    // Progress bar...
    if( (10*iell) / nells != (10*iell+10) / nells )
    {
      printf("Progress Cell integration: %3d%%, ell: %.3e\n",(100*iell+100)/nells,ell);
      if(iell == nells-1) std::cout << std::endl;
    }
    // Set up ODE rhs for this ell
    ODEFunction dC_ell_dlogk = [&] (double logk, const double *y, double *dydx)
    {
      const double k = exp(logk);
      dydx[0] = 4*M_PI*primordial_power_spectrum(k)*f_ell_spline(ell,k)*g_ell_spline(ell,k);
      return GSL_SUCCESS;
    };
    // Solve the ODE
    ODE_C_ell.solve(dC_ell_dlogk,log_k_array,C_ell_IC);
    // Extract end result
    result[iell] = ODE_C_ell.get_final_data_by_component(0);
  }
  Utils::EndTiming("C_ell integration");
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
  fp W15 "ell";
  fp W15 "C_ell_TT";
  fp W15 "C_ell_SW";
  fp W15 "C_ell_ISW";
  fp W15 "C_ell_Doppler";
  fp W15 "C_ell_Quad";
  fp W15 "C_ell_g_tilde";
  fp << "\n";
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    // double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
    //   * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    // double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp W15 ell;
    fp W15 cell_TT_spline(ell) * normfactor;
    fp W15 cell_SW_spline(ell) * normfactor;
    fp W15 cell_ISW_spline(ell) * normfactor;
    fp W15 cell_Doppler_spline(ell) * normfactor;
    fp W15 cell_Quad_spline(ell) * normfactor;
    fp W15 cell_g_tilde_spline(ell) * normfactor;
    // if(Constants.polarization){
    //   fp << cell_EE_spline( ell ) * normfactor  << " ";
    //   fp << cell_TE_spline( ell ) * normfactor  << " ";
    // }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

