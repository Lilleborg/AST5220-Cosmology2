#include"RecombinationHistory.h"
// #include <algorithm>

//====================================================
// Constructors
//====================================================
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================
void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================
void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start,x_end,npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr  = Xe_arr;

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){
    // Get X_e from solving the Saha equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      // Store the result from the Saha equation
      Xe_arr.at(i) = Xe_current;
      ne_arr.at(i) = ne_current;

    } else {
      // Get the rest of Xe from Peebles equation

      // The right hand side of Peebles ODE equation
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      // array of x-values from current to today used in ODEsolver, 
      // init with invalid value (large number) to controll that all elements are overwritten correctly in the following loop
      // there should be no elements with value 100 left in the filled array. (could use linspace here, but this is slightly faster!)
      Vector peebles_x_array(npts_rec_arrays-i,100);
      for (int j = 0; j < npts_rec_arrays-i; j++)
      {
        peebles_x_array.at(j) = x_array.at(j+i);
      }

      ODESolver peebles_Xe_ode;
      Vector peebles_Xe_init{Xe_current};   // Initial condition from the last Xe value found from Saha
      peebles_Xe_ode.solve(dXedx,peebles_x_array,peebles_Xe_init);
      auto Xe_all_data = peebles_Xe_ode.get_data();

      // Fetch results
      for (int j = 0; j < npts_rec_arrays-i; j++)
      {
        Xe_arr.at(j+i) = Xe_all_data[j][0];
        ne_arr.at(j+i) = Xe_arr.at(j+i)*get_number_density_baryons(peebles_x_array.at(j));
      }
      // Solving the full equation until today, so breaking loop
      break;
    }
  }
  
  // Spline the result. Used in get Xe_of_x and ne_of_x methods
  Vector log_Xe_arr = log(Xe_arr);
  Vector log_ne_arr = log(ne_arr);
  log_Xe_of_x_spline.create(x_array,log_Xe_arr);
  log_ne_of_x_spline.create(x_array,log_ne_arr);

  Utils::EndTiming("Xe");
}

//====================================================
// The Saha equation to get ne and Xe in the first regime
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double T_B         = cosmo->get_TCMB()/a;  // Baryon temperature approximation
  const double nb          = get_number_density_baryons(x);

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;

  // Right hand side of Saha equation
  const double temporary_factor = m_e*T_B*k_b/(2*M_PI*hbar*hbar);
  const double F = 1/nb*temporary_factor*sqrt(temporary_factor)*exp(-epsilon_0/k_b/T_B);
  
  // Calculate Xe
  // Determine if we have to use the Taylor approximation in the second order equation
  if (4.0/F < 1e-9){
    Xe = 1.0;
  } else {
    Xe = (-F + F*sqrt(1+4/F))/2;
  }
  
  // Calculate ne
  ne = Xe * nb;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE to get ne and Xe in the first regime
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  
  // Cosmological parameters
  const double T_B         = cosmo->get_TCMB()/a;  // Baryon temperature approximation
  const double nb          = get_number_density_baryons(x);
  const double H           = cosmo->H_of_x(x);
  
  // Some constants to save FLOPS
  const double hbar_square = hbar*hbar;
  const double TB_k_b = k_b*T_B;

  // Expression needed in dXedx
  const double phi_2 = 0.448*log(epsilon_0/TB_k_b);
  const double alpha_2 = 64.0*M_PI/sqrt(27.0*M_PI)*sqrt(epsilon_0/TB_k_b)*phi_2*sigma_T*c*3/(8*M_PI);

  const double beta_factor = m_e*TB_k_b/(2.0*M_PI*hbar_square);
  const double beta = alpha_2*beta_factor*sqrt(beta_factor)*exp(-epsilon_0/TB_k_b);
  const double beta_2 = alpha_2*beta_factor*sqrt(beta_factor)*exp(-epsilon_0/(4.0*TB_k_b));
  // Writing out the full expression for beta into beta_2 to take care of the exponential with possible overflow

  const double n_1s = (1.0-X_e)*nb;
  const double lambda_alpha = H*27.0*epsilon_0*epsilon_0*epsilon_0/(64.0*M_PI*M_PI*n_1s*hbar_square*hbar*c*c*c);
  const double C_r = (lambda_2s1s+lambda_alpha)/(lambda_2s1s+lambda_alpha+beta_2);
  
  dXedx[0] = C_r/H*(beta*(1-X_e) - nb*alpha_2*X_e*X_e);

  // std::cout << "x            " << x << "\n";
  // std::cout << "phi_2        " << phi_2 << "\n";
  // std::cout << "alpha_2      " << alpha_2 << "\n";
  // std::cout << "Current Xe   " << X_e << "\n";
  // std::cout << "beta         " << beta << "\n";
  // std::cout << "beta_2       " << beta_2 << "\n";
  // std::cout << "lambda 2s    " << lambda_2s1s << "\n";
  // std::cout << "C_r          " << C_r << "\n";
  // std::cout << "\n";
  // std::cout << "nb           " << nb << "\n";
  // std::cout << "n_1s         " << n_1s << "\n";
  // std::cout << "lambda alpha " << lambda_alpha << "\n";
  // std::cout << "H            " << H << "\n";

  // std::cout << "dXedx        " << dXedx[0] << "\n";

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1e+4;
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector tau(npts);
  Vector dtaudx_arr(npts);
  Vector dtaudx_arr_analytic(npts);
  Vector ddtauddx_arr(npts);
  Vector vis(npts);
  Vector dvisdx_arr(npts);
  Vector ddvisddx_arr(npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    
    dtaudx[0] = - ne_of_x(x)*Constants.sigma_T*Constants.c/cosmo->H_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  // Finding index of x = 0
  std::vector<double>::iterator x_zero_it_low = std::lower_bound(x_array.begin(),x_array.end(),0);
  int id_x_equal_zero = x_zero_it_low-x_array.begin();
  std::cout << "found x0? " << id_x_equal_zero << " " << x_array[id_x_equal_zero] << "\n";
  // std::vector<double>::iterator x_zero_it_high = std::upper_bound(x_array.begin(),x_array.end(),0);
  // std::cout << "found x0? " << x_zero_it_high-x_array.begin() << " " << x_array[x_zero_it_high-x_array.begin()] << "\n";

  ODESolver tau_ODE;
  Vector tau_init{1e10};
  tau_ODE.solve(dtaudx,x_array,tau_init);
  auto tau_all_data = tau_ODE.get_data();
  for (int i = 0; i < npts; i++)
  {
    tau.at(i) = tau_all_data.at(i)[0] - tau_all_data.at(id_x_equal_zero)[0];
    dtaudx_arr_analytic.at(i) = - ne_of_x(x_array[i])*Constants.sigma_T*Constants.c/cosmo->H_of_x(x_array[i]);
  }
  


  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  tau_of_x_spline.create(x_array,tau);
  //...

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::get_number_density_baryons(double x) const{

  double a = exp(x);
  double nb = cosmo->get_OmegaB(0.0)*cosmo->get_rho_crit(0.0)/(Constants.m_H*a*a*a);

  return nb;
}

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::Xe_of_x(double x) const{

  double log_Xe = log_Xe_of_x_spline(x);
  double res = exp(log_Xe);

  return res;
}

double RecombinationHistory::ne_of_x(double x) const{

  double log_ne = log_ne_of_x_spline(x);
  double res = exp(log_ne);

  return res;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 10000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    // fp << dtaudx_of_x(x)       << " ";
    // fp << ddtauddx_of_x(x)     << " ";
    // fp << g_tilde_of_x(x)      << " ";
    // fp << dgdx_tilde_of_x(x)   << " ";
    // fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

