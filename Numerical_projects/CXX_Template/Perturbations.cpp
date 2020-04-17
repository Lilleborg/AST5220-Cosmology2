#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Vector x_testing = set_up_x_array_resolution();
  // std::cout << x_testing.size() << "\n";
  // std::cout << x_1700 << " " << x_array_full[idx_x1700] << "\n";
  
  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  // compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){
    // Current value of k
    double k = k_array[ik];

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k )
    {
      printf("Progress pert integration: %3d%%, k-value: %.3e , k per Mpc: %.3e\n",(100*ik+100)/n_k,k,k*Constants.Mpc);
      if(ik == n_k-1) std::cout << std::endl;
    }

    //===================================================================
    // Tight couple regime
    //===================================================================
    // Find value to integrate to and set up x-array for tight couple regime
    std::pair<double,int> tight_couple_time_pair = get_tight_coupling_time_and_index(k);
    double x_end_tc = tight_couple_time_pair.first;
    int idx_tc_transition = tight_couple_time_pair.second;
    Vector x_array_tc = Utils::linspace(x_start,x_end_tc,idx_tc_transition);

    // Debugging
    if (x_end_tc > x_1700)
    {
      std::cout << "x_end_tc is larger than x_1700!\n";
      std::cout << "x_1700 " << x_1700 << " " << x_end_tc << std::endl;
    }
    
    // The tight coupling ODE system
    ODESolver ODE_tc_regime;
    auto y_tc_ini = set_ic(x_start, k); // Initial conditions in the tight coupling regime
    ODEFunction dydx_tc = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };
    ODE_tc_regime.solve(dydx_tc,x_array_tc,y_tc_ini);
    auto y_tc_solutions = ODE_tc_regime.get_data();
    auto y_tc_end       = ODE_tc_regime.get_final_data();

    //===================================================================
    // The full system, after tight coupling
    //===================================================================
    // Set up array for after tight coupling
    Vector x_array_after_tc = Utils::linspace(x_end_tc,x_end,n_x-idx_tc_transition);
    printf("number of points in both x arrays: %d, nx: %d\n",int(x_array_tc.size()+x_array_after_tc.size()),n_x);
    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tc_end, x_end_tc, k);

    // The full ODE system
    // ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
    //   return rhs_full_ode(x, k, y, dydx);
    // };

    // Integrate from x_end_tight -> x_end
    // ...
    // ...
    // ...
    // ...
    // ...

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    //===================================================================
    //...
    //...

  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  // const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  // const bool polarization       = Constants.polarization;
  // const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  // double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  // Set the initial conditions in the tight coupling regime
  double fv = 0;
  double Psi = -double(2/3);
  // if (neutrinos)
  // {
  //   fv = cosmo->get_OmegaNu()/(cosmo->get_OmegaNu()+cosmo->get_OmegaR());
  //   Psi = -1/(1.5 + 0.4*fv);
  // }
  // SET: Scalar quantities (Phi, delta, v, ...)
  Phi       = -(1+0.4*fv)*Psi;
  delta_cdm = -1.5*Psi;
  delta_b   = delta_cdm;
  v_cdm     = -Constants.c*k*Psi/(2*cosmo->Hp_of_x(x));
  v_b       = v_cdm;
  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0]  = -0.5*Psi;
  Theta[1]  = -v_cdm/3;
  
  // // SET: Neutrino perturbations (N_ell)
  // if(neutrinos){
  //   // ...
  //   // ...
  // }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================
Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc_end_point, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  // const int n_ell_thetap        = Constants.n_ell_thetap;
  // const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  // const bool polarization       = Constants.polarization;
  // const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  // const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc_end_point[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc_end_point[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc_end_point[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc_end_point[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc_end_point[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc_end_point[Constants.ind_start_theta_tc];
  // const double *Nu_tc           = &y_tc_end_point[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  // double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  // double *Nu              = &y[Constants.ind_start_nu_tc];
  
  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  Phi        = Phi_tc;
  delta_b    = delta_b_tc;
  v_b        = v_b_tc;
  delta_cdm  = delta_cdm_tc;
  v_cdm      = v_cdm_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0]   = Theta_tc[0];
  Theta[1]   = Theta_tc[1];
  const double ck_over_H_p_dtaudx = Constants.c*k/(cosmo->Hp_of_x(x)*rec->dtaudx_of_x(x));
  Theta[2]   = - 20.0/45*ck_over_H_p_dtaudx*Theta[1];
  for (int l = 3; l < n_ell_theta; l++)
  {
    Theta[l] = - l/(2*l+1) * ck_over_H_p_dtaudx * Theta[l-1];
  }
  // SET: Photon polarization perturbations (Theta_p_ell)
  // if(polarization){
  //   // ...
  //   // ...
  // }

  // // SET: Neutrino perturbations (N_ell)
  // if(neutrinos){
  //   // ...
  //   // ...
  // }

  return y;
}

//====================================================
// The time when tight coupling end and the index of that x-value
//====================================================
 std::pair<double,int> Perturbations::get_tight_coupling_time_and_index(const double k) const{
  
  double tau_prime;
  for (int i = 0; i < n_x; i++)
  {
    tau_prime = -rec->dtaudx_of_x(x_array_full[i]);
    if (tau_prime < 10 || tau_prime < 10*k*Constants.c/cosmo->Hp_of_x(x_array_full[i]) || x_array_full[i] > x_1700)
    {
      return std::pair<double,int>(x_array_full[i-1],int(i-1));
    }
  }
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@\nDidn't find tight couple end time!\n"
    "Returning x_1700!\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  return std::pair<double,int>(x_1700,idx_x1700);
}

//====================================================
// After integrating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================
// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  // const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  // const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  // const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  // double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  // Scale factor
  const double a = exp(x);
  const double a_squared = a*a;

  // Other values from earlier milestones and handy quantities
  const double H_p         = cosmo->Hp_of_x(x);
  const double H_p_squared = H_p*H_p;
  const double dH_pdx      = cosmo->dHpdx_of_x(x);
  const double H_0         = cosmo->get_H0();
  const double H_0_squared = H_0*H_0;
  const double OmegaR0     = cosmo->get_OmegaR();
  const double OmegaCDM0   = cosmo->get_OmegaCDM();
  const double OmegaB0     = cosmo->get_OmegaB();
  const double R           = 4*OmegaR0/(3*OmegaB0*a);

  const double dtaudx      = rec->dtaudx_of_x(x);
  const double ddtauddx    = rec->ddtauddx_of_x(x);

  // Constants used in the expressions
  const double ck = Constants.c*k;
  const double ck_squared = ck*ck;
  const double ck_over_H_p = ck/H_p;

  // Theta_2 from initial condition equation
  const double Theta2 = -20*ck_over_H_p*Theta[1]/(45*dtaudx);
  // Psi follows from Phi and Theta2
  const double Psi = - Phi - 12*H_0_squared/(ck_squared*a_squared)*OmegaR0*Theta2;

  // Fill in the expressions for all the derivatives
  // Scalar quantities that are the same as in normal regime
  dPhidx = Psi - ck_squared*Phi/(3*H_p_squared) + H_0_squared/(2*H_p_squared)
    * ( (OmegaCDM0*delta_cdm + OmegaB0*delta_b)/a + 4*OmegaR0*Theta[0]/a_squared);
  ddelta_cdmdx = ck_over_H_p*v_cdm - 3*dPhidx;
  dv_cdmdx     = - v_cdm - ck_over_H_p*Psi;
  ddelta_bdx   = ck_over_H_p*v_b - 3*dPhidx;
  
  // Photon monopole
  dThetadx[0]  = - ck_over_H_p*Theta[1] - dPhidx;

  // Quantities special for tight coupling, q, dv_bdx and dTheta1dx
  const double q_nom = -((1-R)*dtaudx + (1+R)*ddtauddx)*(3*Theta[1]+v_b) - ck_over_H_p*Psi
    + (1-dH_pdx/H_p)*ck_over_H_p*(-Theta[0]+2*Theta2) - ck_over_H_p*dThetadx[0];
  const double q_den = (1+R)*dtaudx + dH_pdx/H_p - 1;
  const double q     = q_nom/q_den;
  dv_bdx             = 1/(1+R)*(-v_b - ck_over_H_p*Psi + R*(q + ck_over_H_p*(-Theta[0] + 2*Theta2) - ck_over_H_p*Psi));
  dThetadx[1]        = (q-dv_bdx)/3.0;

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================
int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  // double Hp = cosmo->Hp_of_x(x);
  // ...

  // Recombination variables
  // ...

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi, delta, v, ...)
  // ...
  // ...
  // ...

  // SET: Photon multipoles (Theta_ell)
  // ...
  // ...
  // ...

  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){
    // ...
    // ...
    // ...
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // ...
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

// Vector Perturbations::set_up_x_array_resolution() const{
//   // Find x-value when recombination is finished following Calin
//   Vector x_array_linspaced = Utils::linspace(x_start,x_end,n_x);
//   bool passed_x_star = false;
//   for (int i = 0; i < n_x; i++)
//   {
    
//   }
  
  
//   Vector x_array_recombination_focus(n_x);
//   x_array_recombination_focus[0] = x_start;
//   int n_x_before = 300;
//   double const dx_before = (x_1700-x_start)/double(n_x_before-1);
//   double dx = dx_before;
  

//   // int nx_test = 100;
//   // Vector x_array_before = Utils::linspace(x_start,x_1700,nx_test/10);
//   // x_array_before.pop_back(); // Removes last entry as will be inserted later
//   // Vector x_array_during = Utils::linspace(x_1700,x_1700*0.9,nx_test/10*8);
//   // x_array_during.pop_back();
//   // Vector x_array_after  = Utils::linspace(x_1700*0.9,x_end,nx_test/10);
//   // x_array = x_array_before;
//   // x_array.insert(x_array.end(), x_array_during.begin(), x_array_during.end());
//   // x_array.insert(x_array.end(), x_array_after.begin(), x_array_after.end());
//   return x_array;
// }

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "x_start:            " << x_start                << "\n";
  std::cout << "x_end:              " << x_end                  << "\n";
  std::cout << "n_x:                " << n_x                    << "\n";
  std::cout << "k_min (1/Mpc):      " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc):      " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:                " << n_k              << "\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";

  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    // fp << get_Theta(x,k,0)   << " ";
    // fp << get_Theta(x,k,1)   << " ";
    // fp << get_Theta(x,k,2)   << " ";
    // fp << get_Phi(x,k)       << " ";
    // fp << get_Psi(x,k)       << " ";
    // fp << get_Pi(x,k)        << " ";
    // fp << get_Source_T(x,k)  << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

