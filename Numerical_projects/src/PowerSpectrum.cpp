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
void PowerSpectrum::solve(bool solve_source_components)
{
  // Concatenating x_array from before, during and after recombination to have higher resolution
  Vector rec_times      = rec->get_time_results();
  double x_rec_start    = rec_times[2];
  double x_rec_end      = rec_times[6];
  Vector x_array_before = Utils::linspace(x_min,x_rec_start,15);
  Vector x_array_during = Utils::linspace(x_rec_start,x_rec_end,172);
  Vector x_array_after  = Utils::linspace(x_rec_end,x_max,15);
  Vector x_array        = x_array_before;
  // As to not have overlapping points, insert x_array_during without start and end points
  // x_array_during is so tightly spaced on a small interval that it should not matter
  x_array.insert(x_array.end(),x_array_during.begin()+1,x_array_during.end()-1);
  x_array.insert(x_array.end(),x_array_after.begin(),x_array_after.end());
  n_x = x_array.size();

  for (int i = 0; i < n_x-1; i++)
  {
    if (x_array[i+1]==x_array[i])
    {
      std::cout << "Overlapping x-points\n" << "at " << i << "\n";
    }
  }
  

  //=========================================================================
  // Make splines for j_ell. 
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // Line of sight integration to get Theta_ell(k)
  //=========================================================================
  line_of_sight_integration(solve_source_components,x_array);

  //=========================================================================
  // Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  //=========================================================================
  
  Vector cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline2D, thetaT_ell_of_k_spline2D);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");

  if (solve_source_components)
  {
    treat_source_components = true;
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
Vector2D PowerSpectrum::line_of_sight_integration_single(std::function<double(double,double)> &source_function,Vector x_array)
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
      printf("Progress LOS integration: %3d%%, k-value per Mpc: %.3e\n",(100*ik+100)/n_k,k*Constants.Mpc);
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
void PowerSpectrum::line_of_sight_integration(bool solve_source_components, const Vector x_array)
{
  std::cout << "Solving with the full temperature source function\n";
  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };
  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(source_function_T,x_array);
  // Create thetaT 2D spline
  thetaT_ell_of_k_spline2D.create(ells,k_array,thetaT_ell_of_k,"thetaT_ell_of_k_spline2D");

  if (solve_source_components)
  {
    std::cout << "Solving with the SW source function\n";
    // Make a function returning the source function
    std::function<double(double,double)> source_function_SW = [&](double x, double k){
      return pert->get_Source_SW(x,k);
    };
    // Do the line of sight integration
    Vector2D thetaSW_ell_of_k = line_of_sight_integration_single(source_function_SW,x_array);
    // Create thetaSW 2D spline
    thetaSW_ell_of_k_spline2D.create(ells,k_array,thetaSW_ell_of_k,"thetaSW_ell_of_k_spline2D");

    std::cout << "Solving with the ISW source function\n";
    // Make a function returning the source function
    std::function<double(double,double)> source_function_ISW = [&](double x, double k){
      return pert->get_Source_ISW(x,k);
    };
    // Do the line of sight integration
    Vector2D thetaISW_ell_of_k = line_of_sight_integration_single(source_function_ISW,x_array);
    // Create thetaISW 2D spline
    thetaISW_ell_of_k_spline2D.create(ells,k_array,thetaISW_ell_of_k,"thetaISW_ell_of_k_spline2D");

    std::cout << "Solving with the Doppler source function\n";
    // Make a function returning the source function
    std::function<double(double,double)> source_function_Doppler = [&](double x, double k){
      return pert->get_Source_Doppler(x,k);
    };
    // Do the line of sight integration
    Vector2D thetaDoppler_ell_of_k = line_of_sight_integration_single(source_function_Doppler,x_array);
    // Create thetaDoppler 2D spline
    thetaDoppler_ell_of_k_spline2D.create(ells,k_array,thetaDoppler_ell_of_k,"thetaDoppler_ell_of_k_spline2D");

    std::cout << "Solving with the Quad source function\n";
    // Make a function returning the source function
    std::function<double(double,double)> source_function_Quad = [&](double x, double k){
      return pert->get_Source_Quad(x,k);
    };
    // Do the line of sight integration
    Vector2D thetaQuad_ell_of_k = line_of_sight_integration_single(source_function_Quad,x_array);
    // Create thetaQuad 2D spline
    thetaQuad_ell_of_k_spline2D.create(ells,k_array,thetaQuad_ell_of_k,"thetaQuad_ell_of_k_spline2D");

    std::cout << "Solving with the g_tilde source function\n";
    // Make a function returning the source function
    std::function<double(double,double)> source_function_g_tilde = [&](double x, double k){
      return rec->g_tilde_of_x(x);
    };
    // Do the line of sight integration
    Vector2D thetag_tilde_ell_of_k = line_of_sight_integration_single(source_function_g_tilde,x_array);
    // Create thetag_tilde 2D spline
    thetag_tilde_ell_of_k_spline2D.create(ells,k_array,thetag_tilde_ell_of_k,"thetag_tilde_ell_of_k_spline2D");
  }
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
double PowerSpectrum::get_component_power_spectrum(const std::string component, const double x, const double k) const{
  
  const double abs_Delta_comp = std::fabs(get_invariant_density(component,x,k));
  const double abs_Delta_comp_squared = abs_Delta_comp*abs_Delta_comp;
  const double pofk = abs_Delta_comp_squared*primordial_power_spectrum(k);

  return pofk;
}

double PowerSpectrum::get_invariant_density(const std::string component, const double x, const double k) const
{
  const double ck = Constants.c*k;
  double delta_comp;
  double omega_comp;
  double v_comp;

  if (component.compare("matter") == 0)
  {
    const double a       = exp(x);
    const double H0      = cosmo->get_H0();
    const double Omega_M = cosmo->get_OmegaCDM()+cosmo->get_OmegaB();
    return ck*ck*pert->get_Phi(x,k)/(1.5*Omega_M/a*H0*H0);;
  }
  else if (component.compare("baryon") == 0)
  {
    delta_comp = pert->get_delta_b(x,k);
    omega_comp = 0;
    v_comp     = pert->get_v_b(x,k);
  }
  else if (component.compare("CDM") == 0)
  {
    delta_comp = pert->get_delta_cdm(x,k);
    omega_comp = 0;
    v_comp     = pert->get_v_cdm(x,k);
  }
  else if (component.compare("radiation") == 0)
  {
    delta_comp = 4.0*pert->get_Theta(x,k,0);
    omega_comp = 1.0/3.0;
    v_comp     = -3*pert->get_Theta(x,k,1);
  }
  else
  {
    std::cout << "Get invariant density didn't find a viable component!!!\nInput was " + component;
  }
  return delta_comp - 3.0*(1+omega_comp)*cosmo->Hp_of_x(x)*v_comp/ck;
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
double PowerSpectrum::get_theta_TT(const double ell, const double k) const
{
  return thetaT_ell_of_k_spline2D(ell,k);
}
//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  std::cout << "Writing output to " << filename << "\n";
  fp W15 "n_k" W15 n_k W15 "n_x" W15 n_x W15 "max ell" W15 ells[ells.size()-1] << "\n";
  fp W15 "ell";
  fp W15 "C_ell_TT";
  if (treat_source_components)
  {
    fp W15 "C_ell_SW";
    fp W15 "C_ell_ISW";
    fp W15 "C_ell_Doppler";
    fp W15 "C_ell_Quad";
    fp W15 "C_ell_g_tilde";
  }
  fp << "\n";
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    fp W15 ell;
    fp W15 cell_TT_spline(ell) * normfactor;
    if (treat_source_components)
    {
      fp W15 cell_SW_spline(ell) * normfactor;
      fp W15 cell_ISW_spline(ell) * normfactor;
      fp W15 cell_Doppler_spline(ell) * normfactor;
      fp W15 cell_Quad_spline(ell) * normfactor;
      fp W15 cell_g_tilde_spline(ell) * normfactor;
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);

}

// Output power spectrum for component
void PowerSpectrum::output_component_power_spectrum(std::vector<std::string> components, std::string filename) const
{
  const double x_eq = cosmo->get_x_equality();
  const double k_eq = exp(x_eq)*cosmo->H_of_x(x_eq)/Constants.c;
  const double Mpc     = Constants.Mpc;
  const double h       = cosmo->get_h();
  const double h_Mpc_3 = h*h*h/Mpc/Mpc/Mpc;
  const double pi_squared_2 = 2*M_PI*M_PI;
  const double prefactor_no_k = h_Mpc_3*pi_squared_2;

  std::ofstream fp(filename.c_str());
  std::cout << "Writing output to " << filename << "\n";
  fp W15 "k_eq" W15 k_eq*Mpc/h << "\n";
  fp W15 "k [h/Mpc],";
  for (int icomp = 0; icomp < components.size(); icomp++)
  {
    fp W15 components[icomp];
    if (icomp != components.size()-1)
    {
      fp << ",";
    }
  }
  fp << "\n";
  
  auto print_data = [&] (const double k)
  {
    const double prefactor = prefactor_no_k/k/k/k;
    fp W15 k*Mpc/h;
    for (int icomp = 0; icomp < components.size(); icomp++)
    {
      fp W15 get_component_power_spectrum(components[icomp],0,k)*prefactor;
    }
    fp << "\n";
  };
  std::for_each(k_array.begin(),k_array.end(),print_data);
}

// Output transfer function and integrand
void PowerSpectrum::output_transfer_integrand(std::string filename, std::vector<int> ell_values) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  std::cout << "Writing output to " << filename << "\n";
  fp W15 "n_k" W15 n_k W15 "n_x" W15 n_x W15 "max ell" W15 ells[ells.size()-1] << "\n";
  fp W15 "ell values:";
  for (int iell = 0; iell < ell_values.size(); iell++)
  {
    fp W15 ell_values[iell];
  }
  fp << "\n";
  
  fp W15 "k";
  fp W15 "theta_1";
  fp W15 "theta_2";
  fp W15 "theta_3";
  fp W15 "theta_4";
  fp W15 "theta_5";
  fp W15 "theta_6";
  fp W15 "theta_1^2/k";
  fp W15 "theta_2^2/k";
  fp W15 "theta_3^2/k";
  fp W15 "theta_4^2/k";
  fp W15 "theta_5^2/k";
  fp W15 "theta_6^2/k";
  fp << "\n";

  auto print_data = [&] (const double k) {
    fp W15 k;
    for (int iell = 0; iell < ell_values.size(); iell++)
    {
      fp W15 get_theta_TT(ell_values[iell],k);
    }
    for (int iell = 0; iell < ell_values.size(); iell++)
    {
      double thetaTT = get_theta_TT(ell_values[iell],k);
      fp W15 thetaTT*thetaTT/k;
    }
    
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);

}