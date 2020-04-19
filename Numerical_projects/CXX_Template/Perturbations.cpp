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

  Vector2D delta_cdm_array_2D (n_x,Vector(n_k));
  Vector2D v_cdm_array_2D     (n_x,Vector(n_k));
  Vector2D delta_b_array_2D   (n_x,Vector(n_k));
  Vector2D v_b_array_2D       (n_x,Vector(n_k));
  Vector2D Phi_array_2D       (n_x,Vector(n_k));
  Vector2D Psi_array_2D       (n_x,Vector(n_k));
  Vector2D Theta0_array_2D    (n_x,Vector(n_k));
  Vector2D Theta1_array_2D    (n_x,Vector(n_k));
  Vector2D Theta2_array_2D    (n_x,Vector(n_k));

  Utils::StartTiming("integrateperturbation");

  // Cosmological parameters
  const double H_0         = cosmo->get_H0();
  const double H_0_squared = H_0*H_0;
  const double OmegaR0     = cosmo->get_OmegaR();

  // Debugging
  bool use_verbose = false;
  
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
    // Solve from x_start -> x_end_tc
    ODE_tc_regime.set_verbose(use_verbose);
    ODE_tc_regime.solve(dydx_tc,x_array_tc,y_tc_ini);
    auto y_tc_solutions = ODE_tc_regime.get_data();
    auto y_tc_end       = ODE_tc_regime.get_final_data();

    //////////////////////////////////////////////////////////////////////////////
    // Store values from tc at start of 2D vectors, using .at() for out of bounds check
    // Constants used in the expressions
    const double ck = Constants.c*k;
    const double ck_squared = ck*ck;

    for (int ix = 0; ix < idx_tc_transition-1; ix++)
    {
      // Scale factor
      const double a = exp(x_array_tc[ix]);
      const double a_squared = a*a;
      const double ck_over_H_p = ck/cosmo->Hp_of_x(x_array_tc[ix]);
      const double dtaudx = rec->dtaudx_of_x(x_array_tc[ix]);

      // Scalar quantities
      delta_cdm_array_2D.at(ix).at(ik) = y_tc_solutions.at(ix).at(Constants.ind_deltacdm_tc);
      v_cdm_array_2D.at(ix).at(ik)     = y_tc_solutions.at(ix).at(Constants.ind_vcdm_tc);
      delta_b_array_2D.at(ix).at(ik)   = y_tc_solutions.at(ix).at(Constants.ind_deltab_tc);
      v_b_array_2D.at(ix).at(ik)       = y_tc_solutions.at(ix).at(Constants.ind_vb_tc);
      Phi_array_2D.at(ix).at(ik)       = y_tc_solutions.at(ix).at(Constants.ind_Phi_tc);
      
      // Multipoles
      Theta0_array_2D.at(ix).at(ik)    = y_tc_solutions.at(ix).at(Constants.ind_start_theta_tc);
      Theta1_array_2D.at(ix).at(ik)    = y_tc_solutions.at(ix).at(Constants.ind_start_theta_tc+1);
      Theta2_array_2D.at(ix).at(ik)    = -20.0*ck_over_H_p/(45.0*dtaudx)*Theta1_array_2D.at(ix).at(ik);
      // Psi not dynamical so not in ODE solution.
      Psi_array_2D.at(ix).at(ik)       = - Phi_array_2D.at(ix).at(ik)
        - 12*H_0_squared/(ck_squared*a_squared)*OmegaR0*Theta2_array_2D.at(ix).at(ik);
    }
    // End storing tight couple data
    //////////////////////////////////////////////////////////////////////////////
    //===================================================================
    // The full system, after tight coupling
    //===================================================================
    // Set up array for after tight coupling
    Vector x_array_after_tc = Utils::linspace(x_end_tc,x_end,n_x-idx_tc_transition);
    if (x_array_tc.back() == x_array_after_tc[0])
    {
      printf("Duplicate x-value before and after tc, %e, index %d\n",x_array_tc.back(),x_array_tc.size());
    }

    if (x_array_tc.size()+x_array_after_tc.size() != n_x)
    {
      printf("number of points in x arrays before and after tc not equal max points: %d, nx: %d\n",int(x_array_tc.size()+x_array_after_tc.size()),n_x);
    }
    // The full ODE system after tight couping
    // Set up initial conditions (y_tc_end is the solution at the end of tight coupling)
    ODESolver ODE_after_tc;
    auto y_after_tc_ini = set_ic_after_tight_coupling(y_tc_end, x_end_tc, k);
    ODEFunction dydx_after_tc = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };
    // Solve from x_end_tight -> x_end
    ODE_after_tc.set_verbose(use_verbose);
    ODE_after_tc.solve(dydx_after_tc,x_array_after_tc,y_after_tc_ini);
    auto y_after_tc_solutions = ODE_after_tc.get_data();

    //////////////////////////////////////////////////////////////////////////////
    // Store values from from after tc at end of 2D vectors, using .at() for out of bounds check
    for (int ix = 0; ix < n_x-(idx_tc_transition); ix++)
    {
      // Scale factor
      const double a = exp(x_array_after_tc[ix]);
      const double a_squared = a*a;
      const double ck_over_H_p = ck/cosmo->Hp_of_x(x_array_after_tc[ix]);
      const double dtaudx = rec->dtaudx_of_x(x_array_after_tc[ix]);

      // Scalar quantities
      delta_cdm_array_2D.at(ix+idx_tc_transition).at(ik) = y_after_tc_solutions.at(ix).at(Constants.ind_deltacdm);
      v_cdm_array_2D.at(ix+idx_tc_transition).at(ik)     = y_after_tc_solutions.at(ix).at(Constants.ind_vcdm);
      delta_b_array_2D.at(ix+idx_tc_transition).at(ik)   = y_after_tc_solutions.at(ix).at(Constants.ind_deltab);
      v_b_array_2D.at(ix+idx_tc_transition).at(ik)       = y_after_tc_solutions.at(ix).at(Constants.ind_vb);
      Phi_array_2D.at(ix+idx_tc_transition).at(ik)       = y_after_tc_solutions.at(ix).at(Constants.ind_Phi);
      
      // Multipoles
      Theta0_array_2D.at(ix+idx_tc_transition).at(ik)    = y_after_tc_solutions.at(ix).at(Constants.ind_start_theta);
      Theta1_array_2D.at(ix+idx_tc_transition).at(ik)    = y_after_tc_solutions.at(ix).at(Constants.ind_start_theta+1);
      Theta2_array_2D.at(ix+idx_tc_transition).at(ik)    = y_after_tc_solutions.at(ix).at(Constants.ind_start_theta+2);
      // Psi not dynamical so not in ODE solution.
      Psi_array_2D.at(ix+idx_tc_transition).at(ik)       = - Phi_array_2D.at(ix+idx_tc_transition).at(ik)
        - 12.0*H_0_squared/(ck_squared*a_squared)*OmegaR0*Theta2_array_2D.at(ix+idx_tc_transition).at(ik);
    }
    // End storing tight couple data
    //////////////////////////////////////////////////////////////////////////////

  }
  Utils::EndTiming("integrateperturbation");

  // Check for duplicate values
  for (int ik = 0; ik < n_k; ik++)
  {
    for (int ix = 0; ix < n_x-1; ix++)
    {
      if (Theta0_array_2D.at(ix).at(ik) == Theta0_array_2D.at(ix+1).at(ik))
      {
        printf("Duplciate value in 2D array, ix: %d, Theta0(ix): %e, Theta0(ix+1): %e\n"\
          ,ix,Theta0_array_2D.at(ix).at(ik),Theta0_array_2D.at(ix+1).at(ik));
      }
    }
  }
  

  //=============================================================================
  // Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  delta_cdm_spline.create(x_array_full,k_array,delta_cdm_array_2D,"delta_cdm_spline");
  v_cdm_spline.create(x_array_full,k_array,v_cdm_array_2D,"v_cdm_spline");
  delta_b_spline.create(x_array_full,k_array,delta_b_array_2D,"delta_b_spline");
  v_b_spline.create(x_array_full,k_array,v_b_array_2D,"v_b_spline");
  Phi_spline.create(x_array_full,k_array,Phi_array_2D,"Phi_spline");
  Psi_spline.create(x_array_full,k_array,Psi_array_2D,"Psi_spline");
  Pi_spline.create(x_array_full,k_array,Theta2_array_2D,"Pi_spline"); // Not including polarization
  
  Theta0_spline.create(x_array_full,k_array,Theta0_array_2D,"Theta0_spline");
  Theta1_spline.create(x_array_full,k_array,Theta1_array_2D,"Theta1_spline");
  Theta2_spline.create(x_array_full,k_array,Theta2_array_2D,"Theta2_spline");
  vector_of_Theta_splines.push_back(Theta0_spline);
  vector_of_Theta_splines.push_back(Theta1_spline);
  vector_of_Theta_splines.push_back(Theta2_spline);
}

//====================================================
// Set IC at the start of the run (this is in the tight coupling regime)
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
  double fv = 0.0;
  double Psi = -2.0/3.0;
  // if (neutrinos)
  // {
  //   fv = cosmo->get_OmegaNu()/(cosmo->get_OmegaNu()+cosmo->get_OmegaR());
  //   Psi = -1/(1.5 + 0.4*fv);
  // }
  // Scalar quantities (Phi, delta, v, ...)
  Phi       = -(1.0+0.4*fv)*Psi;
  delta_cdm = -1.5*Psi;
  delta_b   = delta_cdm;
  v_cdm     = -Constants.c*k*Psi/(2.0*cosmo->Hp_of_x(x));
  v_b       = v_cdm;
  // Photon temperature perturbations (Theta_ell)
  Theta[0]  = -0.5*Psi;
  Theta[1]  = -v_cdm/3.0;
  
  // // SET: Neutrino perturbations (N_ell)
  // if(neutrinos){
  //   // ...
  //   // ...
  // }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling regime ends
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
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];
  // double *Theta_p         = &y[Constants.ind_start_thetap];
  // double *Nu              = &y[Constants.ind_start_nu];
  
  // Scalar quantities (Gravitational potential, baryons and CDM)
  Phi        = Phi_tc;
  delta_b    = delta_b_tc;
  v_b        = v_b_tc;
  delta_cdm  = delta_cdm_tc;
  v_cdm      = v_cdm_tc;

  // Photon temperature perturbations (Theta_ell)
  const double ck_over_H_p_dtaudx = Constants.c*k/(cosmo->Hp_of_x(x)*rec->dtaudx_of_x(x));
  // double l_as_double;
  Theta[0]   = Theta_tc[0];
  Theta[1]   = Theta_tc[1];
  Theta[2]   = - 20.0/45.0*ck_over_H_p_dtaudx*Theta[1];
  for (int l = 3; l < n_ell_theta; l++)
  {
    // l_as_double = double(l);
    Theta[l] = - l/(2.0*l+1.0) * ck_over_H_p_dtaudx * Theta[l-1];
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
// The right hand side of the perturbations ODE in the tight coupling regime
//====================================================
// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){
  
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

  // Cosmological parameters and variables
  const double H_p         = cosmo->Hp_of_x(x);
  const double H_p_squared = H_p*H_p;
  const double dH_pdx      = cosmo->dHpdx_of_x(x);
  const double H_0         = cosmo->get_H0();
  const double H_0_squared = H_0*H_0;
  const double OmegaR0     = cosmo->get_OmegaR();
  const double OmegaCDM0   = cosmo->get_OmegaCDM();
  const double OmegaB0     = cosmo->get_OmegaB();
  const double R           = 4.0*OmegaR0/(3.0*OmegaB0*a);

  // Recombination variables
  const double dtaudx      = rec->dtaudx_of_x(x);
  const double ddtauddx    = rec->ddtauddx_of_x(x);

  // Constants used in the expressions
  const double ck = Constants.c*k;
  const double ck_squared = ck*ck;
  const double ck_over_H_p = ck/H_p;

  // Theta_2 from initial condition equation
  const double Theta2 = -20.0*ck_over_H_p/(45.0*dtaudx)*Theta[1];
  // Psi follows from Phi and Theta2
  const double Psi = - Phi - 12.0*H_0_squared/(ck_squared*a_squared)*OmegaR0*Theta2;

  // Fill in the expressions for all the derivatives
  // Scalar quantities that are the same as in normal regime
  dPhidx       = Psi - ck_squared*Phi/(3.0*H_p_squared) + H_0_squared/(2.0*H_p_squared)
    * ( (OmegaCDM0*delta_cdm + OmegaB0*delta_b)/a + 4.0*OmegaR0*Theta[0]/a_squared);
  ddelta_cdmdx = ck_over_H_p*v_cdm - 3.0*dPhidx;
  dv_cdmdx     = - v_cdm - ck_over_H_p*Psi;
  ddelta_bdx   = ck_over_H_p*v_b - 3.0*dPhidx;
  
  // Photon monopole
  dThetadx[0]  = - ck_over_H_p*Theta[1] - dPhidx;

  // Quantities special for tight coupling, q, dv_bdx and dTheta1dx
  const double q_nom = -((1.0-R)*dtaudx + (1.0+R)*ddtauddx)*(3.0*Theta[1]+v_b) - ck_over_H_p*Psi
    + (1.0-dH_pdx/H_p)*ck_over_H_p*(-Theta[0]+2.0*Theta2) - ck_over_H_p*dThetadx[0];
  const double q_den = (1.0+R)*dtaudx + dH_pdx/H_p - 1.0;
  const double q     = q_nom/q_den;
  dv_bdx             = 1.0/(1.0+R)*(-v_b - ck_over_H_p*Psi + R*(q + ck_over_H_p*(-Theta[0] + 2.0*Theta2) - ck_over_H_p*Psi));
  dThetadx[1]        = (q-dv_bdx)/3.0;

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================
int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  // const int n_ell_thetap        = Constants.n_ell_thetap;
  // const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  // const bool polarization       = Constants.polarization;
  // const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  // const double *Theta_p         = &y[Constants.ind_start_thetap];
  // const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  // double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  // double *dNudx           = &dydx[Constants.ind_start_nu];

  // Scale factor
  const double a = exp(x);
  const double a_squared = a*a;

  // Cosmological parameters and variables
  const double H_p         = cosmo->Hp_of_x(x);
  const double H_p_squared = H_p*H_p;
  // const double dH_pdx      = cosmo->dHpdx_of_x(x);
  const double H_0         = cosmo->get_H0();
  const double H_0_squared = H_0*H_0;
  const double OmegaR0     = cosmo->get_OmegaR();
  const double OmegaCDM0   = cosmo->get_OmegaCDM();
  const double OmegaB0     = cosmo->get_OmegaB();
  const double R           = 4.0*OmegaR0/(3.0*OmegaB0*a);

  // Recombination variables
  const double dtaudx      = rec->dtaudx_of_x(x);

  // Constants used in the expressions
  const double ck = Constants.c*k;
  const double ck_squared = ck*ck;
  const double ck_over_H_p = ck/H_p;

  // Psi follows from Phi and Theta2
  const double Psi = - Phi - 12.0*H_0_squared/(ck_squared*a_squared)*OmegaR0*Theta[2];

  // Scalar quantities (Phi, delta, v, ...)
  dPhidx       = Psi - ck_squared*Phi/(3.0*H_p_squared) + H_0_squared/(2.0*H_p_squared)
    * ( (OmegaCDM0*delta_cdm + OmegaB0*delta_b)/a + 4.0*OmegaR0*Theta[0]/a_squared);
  ddelta_cdmdx = ck_over_H_p*v_cdm - 3.0*dPhidx;
  dv_cdmdx     = - v_cdm - ck_over_H_p*Psi;
  ddelta_bdx   = ck_over_H_p*v_b - 3.0*dPhidx;
  dv_bdx       = - v_b - ck_over_H_p*Psi + dtaudx*R*(3.0*Theta[1]+v_b); // Differ from tc

  // Photon multipoles (Theta_ell)
  dThetadx[0]  = - ck_over_H_p*Theta[1] - dPhidx;
  dThetadx[1]  = ck_over_H_p/3.0*(Theta[0]+Psi) - 2.0*ck_over_H_p/3.0*Theta[2] + dtaudx*(Theta[1]+v_b/3.0); // Differ from tc
  for (int l = 2; l < n_ell_theta-1; l++) // loop from l=2 to l<l_max (l_max=n_ell_theta-1)
  {
    dThetadx[l] = l*ck_over_H_p/(2.0*l+1.0)*Theta[l-1] - (l+1.0)*ck_over_H_p/(2.0*l+1.0)*Theta[l+1]
      + dtaudx*Theta[l];
    if (l==2) // Add value if l==2, due to delta function in expression
    {
      dThetadx[l] -= dtaudx*Theta[2]/10.0;
    }
  }
  dThetadx[n_ell_theta] = ck_over_H_p*Theta[n_ell_theta-1]
  + Theta[n_ell_theta]*(dtaudx - Constants.c*(n_ell_theta+1.0)/(H_p*cosmo->eta_of_x(x)));


  // SET: Photon polarization multipoles (Theta_p_ell)
  // if(polarization){
  //   // ...
  //   // ...
  //   // ...
  // }

  // // SET: Neutrino mutlipoles (Nu_ell)
  // if(neutrinos){
  //   // ...
  //   // ...
  //   // ...
  // }

  return GSL_SUCCESS;
}

//====================================================
// After integrating the perturbation compute the source function(s)
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
  // if(Constants.polarization){
  //   SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  // }

  Utils::EndTiming("source");
}

//====================================================
// Get methods
//====================================================

//====================================================
// The time when tight coupling end and the index of that x-value
//====================================================
std::pair<double,int> Perturbations::get_tight_coupling_time_and_index(const double k) const{
  
  double tau_prime;
  for (int i = 0; i < n_x; i++)
  {
    // dtaudx always negative, instead of using absolute value just force it to be positive
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
// double Perturbations::get_Source_E(const double x, const double k) const{
//   return SE_spline(x,k);
// }
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return vector_of_Theta_splines.at(ell)(x,k);
}
// double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
//   return Theta_p_spline[ell](x,k);
// }
// double Perturbations::get_Nu(const double x, const double k, const int ell) const{
//   return Nu_spline[ell](x,k);
// }

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
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    // fp << get_Source_T(x,k)  << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

