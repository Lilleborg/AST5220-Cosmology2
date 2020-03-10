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
  
  // Set up x-array and make arrays to store X_e(x) and n_e(x) in
  Vector x_array = Utils::linspace(x_start,x_end,npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_arr_only_Saha = Xe_arr;
  Vector ne_arr  = Xe_arr;

  // Save index of when leaving Saha regime and start using Peebles
  int idx_Peebles_transition = 0;
  double Xe_Peebles_transition = 0;

  // Calculate recombination history using Saha approximation all the way
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++)
  {
    // Get X_e and ne from solving the Saha equation
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // If we enter the Peebles regime for the first time, keep track of index i and value of Xe_current
    if(Xe_current < Xe_saha_limit && saha_regime == true)
    {
      saha_regime = false;
      idx_Peebles_transition = i;
      Xe_Peebles_transition = Xe_current;
    }

    // Store the result from the Saha equation
    Xe_arr[i] = Xe_current; // Keep all solutions, overwrite the ones who are in Peebles regime later
    ne_arr[i] = ne_current;
    Xe_arr_only_Saha[i] = Xe_current; // Keep the solutions obtained by Saha eq only
  }
  
  ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
    return rhs_peebles_ode(x, Xe, dXedx);
  };
  
  // array of x-values from current to today used in ODEsolver, 
  // init with invalid value (large number) to controll that all elements are overwritten correctly in the following loop
  // there should be no elements with value 100 left in the filled array.
  // (could have used linspace here, but this way is slightly faster!)
  int i = idx_Peebles_transition;
  Vector peebles_x_array(npts_rec_arrays-i,100);
  for (int j = 0; j < npts_rec_arrays-i; j++)
  {
    peebles_x_array[j] = x_array[j+i];
  }

  ODESolver peebles_Xe_ode;
  Vector peebles_Xe_init{Xe_Peebles_transition};   // Initial condition from the last Xe value found from Saha
  peebles_Xe_ode.solve(dXedx,peebles_x_array,peebles_Xe_init);
  auto peebles_Xe_solution = peebles_Xe_ode.get_data_by_component(0);

  // Store results
  for (int j = 0; j < npts_rec_arrays-i; j++)
  {
    // Overwrite the solutions found using Saha, with the proper Peebles solution
    Xe_arr[j+i] = peebles_Xe_solution[j];
    ne_arr[j+i] = Xe_arr[j+i]*get_number_density_baryons(peebles_x_array[j]);
  }

  // Spline the result. Used in get Xe_of_x and ne_of_x methods
  Vector log_Xe_arr = log(Xe_arr);
  Vector log_Xe_arr_only_Saha = log(Xe_arr_only_Saha);
  Vector log_ne_arr = log(ne_arr);
  log_Xe_of_x_spline.create(x_array,log_Xe_arr,"log Xe");
  log_Xe_of_x_spline_only_Saha.create(x_array,log_Xe_arr_only_Saha,"log Xe Saha");
  log_ne_of_x_spline.create(x_array,log_ne_arr,"log ne");
}

//====================================================
// The Saha equation to get ne and Xe in the first regime
//====================================================
Doublepair RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
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
  const double F = 1/nb*temporary_factor*sqrt(temporary_factor)*exp(-epsilon_0/(k_b*T_B));
  
  // Calculate Xe
  // If near endpoint, take care of instability and set solution to basically zero (splines the log so cant put it to zero)
  if (F < 1e-20)
    Xe = 1e-20;
  // Determine if we have to use the Taylor approximation in the second order equation
  else if (F > 1e+9)
    Xe = 1.0;
  else
    Xe = F*(-1 + sqrt(1.0+4.0/F))/2.0;
  
  // Calculate ne
  ne = Xe * nb;

  return Doublepair(Xe, ne);
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
  // const double hbar_cube = hbar_square*hbar;
  // const double ccc = c*c*c;
  const double TB_k_b = k_b*T_B;

  // Expression needed in dXedx
  const double phi_2 = 0.448*log(epsilon_0/TB_k_b);
  const double alpha_2 = 8/sqrt(3.0*M_PI)*sqrt(epsilon_0/TB_k_b)*phi_2*sigma_T*c;

  const double beta_factor = m_e*TB_k_b/(2.0*M_PI*hbar_square);
  const double beta = alpha_2*beta_factor*sqrt(beta_factor)*exp(-epsilon_0/TB_k_b);
  const double beta_2 = alpha_2*beta_factor*sqrt(beta_factor)*exp(-epsilon_0/(4.0*TB_k_b));
  // Writing out the full expression for beta into beta_2 to take care of the exponential with possible overflow

  const double n_1s = (1.0-X_e)*nb;
  const double lambda_alpha = H*27.0*epsilon_0*epsilon_0*epsilon_0/(64.0*M_PI*M_PI*n_1s*hbar_square*hbar*c*c*c);
  const double C_r = (lambda_2s1s+lambda_alpha)/(lambda_2s1s+lambda_alpha+beta_2);
  
  dXedx[0] = C_r/H*(beta*(1-X_e) - nb*alpha_2*X_e*X_e);

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================
void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over and initiate arrays to store values
  Vector x_array = Utils::linspace(x_start, x_end, npts_tau);
  Vector tau_arr(npts_tau);
  Vector dtaudx_arr(npts_tau);
  Vector vis(npts_tau);
  Vector dvisdx_arr(npts_tau);

  // Set up and solve the ODE for tau
  // The ODE right hand side of dtaudx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    dtaudx[0] = - ne_of_x(x)*Constants.sigma_T*Constants.c/cosmo->H_of_x(x);
    return GSL_SUCCESS;
  };

  // Finding index of x = 0
  std::vector<double>::iterator x_zero_it_low = std::lower_bound(x_array.begin(),x_array.end(),0);
  int id_x_equal_zero = x_zero_it_low-x_array.begin();

  ODESolver tau_ODE;
  Vector tau_init{1e3}; // Some initial value, normalzed later in loop as we know tau(x=0) = 0
  tau_ODE.solve(dtaudx,x_array,tau_init,gsl_odeiv2_step_rkf45);
  // Fetch results
  Vector tau_data = tau_ODE.get_data_by_component(0); // Not normalized yet!
  dtaudx_arr = tau_ODE.get_derivative_data_by_component(0);
  
  // Spline the derivative array, used in the first derivative of the visibility function
  tau_deriv_of_x_spline.create(x_array,dtaudx_arr,"tau derivative");
  
  for (int i = 0; i < npts_tau; i++)
  {
    tau_arr[i] = tau_data[i] - tau_data[id_x_equal_zero]; // normalise with today value
    vis[i] = - dtaudx_arr[i]*exp(-tau_arr[i]);                  // store visibility function
    // Store its first derivative using the product and chain rule
    // This way the second derivative of visibility func can be obtained from the first derivative spline
    dvisdx_arr[i] = exp(-tau_arr[i])*(dtaudx_arr[i]*dtaudx_arr[i]-ddtauddx_of_x(x_array[i]));
  }
  
  // Spline the tau results
  tau_of_x_spline.create(x_array,tau_arr,"tau");

  // Spline the visibility resutls
  g_tilde_of_x_spline.create(x_array,vis,"g tilde");
  g_tilde_deriv_of_x_spline.create(x_array,dvisdx_arr,"g tilde deriv");

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================
// Get times for last scattering, in x and z
Vector RecombinationHistory::get_time_results() const{
  Vector res(6);
  // Using the tau spline and binary search for value method to find tau = 1
  Doublepair xrange(-10.0,-5.0);  // Range of x-value to search in
  // x and z value when tau equals one stored
  res[0] = Utils::binary_search_for_value(tau_of_x_spline,1.0,xrange);
  res[1] = 1/exp(res[0]) - 1;

  // Using Xe spline to search for Xe = 0.5, the spline is log so search for log(0.5)
  res[2] = Utils::binary_search_for_value(log_Xe_of_x_spline,log(0.5),xrange);
  res[3] = 1/exp(res[2]) - 1;

  // Using Xe_saha_only spline to search for Xe = 0.5, the spline is log so search for log(0.5)
  res[4] = Utils::binary_search_for_value(log_Xe_of_x_spline_only_Saha,log(0.5),xrange);
  res[5] = 1/exp(res[2]) - 1;
  return res;
}

// Get number density of baryons from cosmo object
double RecombinationHistory::get_number_density_baryons(double x) const{
  double a = exp(x);
  double nb = cosmo->get_OmegaB(0.0)*cosmo->get_rho_crit(0.0)/(Constants.m_H*a*a*a);

  return nb;
}

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

// Returns the derivative of tau using the tau_deriv_spline
double RecombinationHistory::dtaudx_of_x(double x) const{
  return tau_deriv_of_x_spline(x);
}

// Returns the second derivative of tau using the tau_deriv_spline.deriv_x
double RecombinationHistory::ddtauddx_of_x(double x) const{
  return tau_deriv_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

// Returns the derivative of g using the g_tilde_deriv_spline.deriv_x
double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  // return g_tilde_of_x_spline.deriv_x(x);
  return g_tilde_deriv_of_x_spline(x);
}

// Returns the second derivative of g using the g_tilde_deriv_spline.deriv_xx
double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  // return g_tilde_of_x_spline.deriv_xx(x);
  return g_tilde_deriv_of_x_spline.deriv_x(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  double log_Xe = log_Xe_of_x_spline(x);
  double res = exp(log_Xe);

  return res;
}

double RecombinationHistory::Xe_of_x_Saha_approx(double x) const{
  double log_Xe = log_Xe_of_x_spline_only_Saha(x);
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

// Print values of x and z for last scattering, using tau(x) = tau(z) = 1;
void RecombinationHistory::print_time_results() const{
  Vector times = get_time_results();

  std::cout << "\nTime for last scattering,\ntau(x_star) = tau(z_star) = 1:\n";
  std::cout << "x_star:      " << times[0] << "\n";
  std::cout << "z_star:      " << times[1] << "\n";
  
  std::cout << "\nTime for half-way recombination,\nXe(x_rec) = Xe(z_rec) = 0.5:\n";
  std::cout << "x_rec:      " << times[2] << "\n";
  std::cout << "z_rec:      " << times[3] << "\n";
  
  std::cout << "\nTime for half-way recombination using only Saha approximation:\n";
  std::cout << "x_rec_Saha:      " << times[4] << "\n";
  std::cout << "z_rec_Saha:      " << times[5] << "\n";
  std::cout << std::endl;
}

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  // Create x_array to write to file using the splines. Choose a narrower interval than 
  // the one used to solve the equations to avoid boundary problems
  const int npts       = 1e+4;
  const double x_min   = x_start+3;
  const double x_max   = x_end-2;
  Vector x_array = Utils::linspace(x_min, x_max, npts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                      << " ";
    fp << Xe_of_x(x)             << " ";
    fp << Xe_of_x_Saha_approx(x) << " ";
    fp << ne_of_x(x)             << " ";
    fp << tau_of_x(x)            << " ";
    fp << dtaudx_of_x(x)         << " ";
    fp << ddtauddx_of_x(x)       << " ";
    fp << g_tilde_of_x(x)        << " ";
    fp << dgdx_tilde_of_x(x)     << " ";
    fp << ddgddx_tilde_of_x(x)   << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

