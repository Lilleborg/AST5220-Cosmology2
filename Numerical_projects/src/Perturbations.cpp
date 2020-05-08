#include"Perturbations.h"


#include <sstream>
#include <iomanip>


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
  
  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================
void Perturbations::integrate_perturbations(){

  // Scalar quantities
  Vector2D delta_cdm_array_2D  (n_x,Vector(n_k));
  Vector2D v_cdm_array_2D      (n_x,Vector(n_k));
  Vector2D delta_b_array_2D    (n_x,Vector(n_k));
  Vector2D v_b_array_2D        (n_x,Vector(n_k));
  Vector2D dv_b_dx_array_2D    (n_x,Vector(n_k));
  Vector2D Phi_array_2D        (n_x,Vector(n_k));
  Vector2D dPhi_dx_array_2D    (n_x,Vector(n_k));
  Vector2D Psi_array_2D        (n_x,Vector(n_k));
  // Multipoles and their derivatives
  Vector2D Theta0_array_2D     (n_x,Vector(n_k));
  Vector2D dTheta0_dx_array_2D (n_x,Vector(n_k));
  Vector2D Theta1_array_2D     (n_x,Vector(n_k));
  Vector2D dTheta1_dx_array_2D (n_x,Vector(n_k));
  Vector2D Theta2_array_2D     (n_x,Vector(n_k));
  Vector2D Theta3_array_2D     (n_x,Vector(n_k));

  // Need some derivative quantities in calculating the source function
  // Start simple, use spline to get dTheta2dx and dTheta3dx as no expression given for this in TC
  // might not need dtheta3dx if using approximation Pi = Theta2 in next milestone

  // Cosmological parameters
  const double H_0         = cosmo->get_H0();
  const double H_0_squared = H_0*H_0;
  const double OmegaR0     = cosmo->get_OmegaR();

  Utils::StartTiming("integrateperturbation");
  
  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){
    // Current value of k
    double k = k_array[ik];

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k )
    {
      printf("Progress pert integration: %3d%%, k-value: %.3e , k per Mpc: %.5f\n",(100*ik+100)/n_k,k,k*Constants.Mpc);
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
    
    // The tight coupling ODE system
    ODESolver ODE_tc_regime;
    auto y_tc_ini = set_ic(x_start, k); // Initial conditions in the tight coupling regime
    
    ODEFunction dydx_tc = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };
    // Solve from x_start -> x_end_tc
    ODE_tc_regime.solve(dydx_tc,x_array_tc,y_tc_ini);
    auto y_tc_solutions = ODE_tc_regime.get_data();
    auto y_deriv_tc_sol = ODE_tc_regime.get_derivative_data();
    auto y_tc_end       = ODE_tc_regime.get_final_data();

    //===================================================================
    // The full system, after tight coupling
    //===================================================================
    // Set up array for after tight coupling
    Vector x_array_after_tc = Utils::linspace(x_end_tc,x_end,n_x-idx_tc_transition+1);

    // The full ODE system after tight couping
    // Set up initial conditions (y_tc_end is the solution at the end of tight coupling)
    ODESolver ODE_after_tc;
    auto y_after_tc_ini = set_ic_after_tight_coupling(y_tc_end, x_end_tc, k);
    ODEFunction dydx_after_tc = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };
    // Solve from x_end_tight -> x_end
    ODE_after_tc.solve(dydx_after_tc,x_array_after_tc,y_after_tc_ini);
    auto y_after_tc_solutions = ODE_after_tc.get_data();
    auto y_deriv_after_tc_sol = ODE_after_tc.get_derivative_data();

    //////////////////////////////////////////////////////////////////////////////
    // Store values from tc at start of 2D vectors
    // Constants used in the expressions
    const double ck = Constants.c*k;
    const double ck_squared = ck*ck;

    for (int ix = 0; ix < idx_tc_transition; ix++)
    {
      // Scale factor
      const double a = exp(x_array_tc[ix]);
      const double a_squared = a*a;
      const double ck_over_H_p = ck/cosmo->Hp_of_x(x_array_tc[ix]);
      const double dtaudx = rec->dtaudx_of_x(x_array_tc[ix]);

      // Scalar quantities
      delta_cdm_array_2D[ix][ik] = y_tc_solutions[ix][Constants.ind_deltacdm_tc];
      v_cdm_array_2D[ix][ik]     = y_tc_solutions[ix][Constants.ind_vcdm_tc];
      delta_b_array_2D[ix][ik]   = y_tc_solutions[ix][Constants.ind_deltab_tc];
      v_b_array_2D[ix][ik]       = y_tc_solutions[ix][Constants.ind_vb_tc];
      Phi_array_2D[ix][ik]       = y_tc_solutions[ix][Constants.ind_Phi_tc];
      
      // Multipoles
      Theta0_array_2D[ix][ik]    = y_tc_solutions[ix][Constants.ind_start_theta_tc];
      Theta1_array_2D[ix][ik]    = y_tc_solutions[ix][Constants.ind_start_theta_tc+1];
      Theta2_array_2D[ix][ik]    = -20.0*ck_over_H_p/(45.0*dtaudx)*Theta1_array_2D[ix][ik];
      // Using initial condition with l = 3:
      Theta3_array_2D[ix][ik]    = -3.0/(7.0)*ck_over_H_p/dtaudx*Theta2_array_2D[ix][ik];
      // Derivatives
      dTheta0_dx_array_2D[ix][ik]= y_deriv_tc_sol[ix][Constants.ind_start_theta_tc];
      dTheta1_dx_array_2D[ix][ik]= y_deriv_tc_sol[ix][Constants.ind_start_theta_tc+1];
      dPhi_dx_array_2D[ix][ik]   = y_deriv_tc_sol[ix][Constants.ind_Phi_tc];
      dv_b_dx_array_2D[ix][ik]   = y_deriv_tc_sol[ix][Constants.ind_vb_tc];

      // Maybe implement expressions for dTheta2dx and dTheta3dx if have time in next milestone
      // for now get them from splines

      // Psi not dynamical so not in ODE solution.
      Psi_array_2D[ix][ik]       = - Phi_array_2D[ix][ik]
        - 12*H_0_squared/(ck_squared*a_squared)*OmegaR0*Theta2_array_2D[ix][ik];
    }
    // End storing tight couple data
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // Store values from after tc at end of 2D vectors
    for (int ix = 1; ix < n_x-(idx_tc_transition)+1; ix++)
    {
      // Scale factor
      const double a = exp(x_array_after_tc[ix]);
      const double a_squared = a*a;

      // Scalar quantities
      delta_cdm_array_2D[ix-1+idx_tc_transition][ik] = y_after_tc_solutions[ix][Constants.ind_deltacdm];
      v_cdm_array_2D[ix-1+idx_tc_transition][ik]     = y_after_tc_solutions[ix][Constants.ind_vcdm];
      delta_b_array_2D[ix-1+idx_tc_transition][ik]   = y_after_tc_solutions[ix][Constants.ind_deltab];
      v_b_array_2D[ix-1+idx_tc_transition][ik]       = y_after_tc_solutions[ix][Constants.ind_vb];
      Phi_array_2D[ix-1+idx_tc_transition][ik]       = y_after_tc_solutions[ix][Constants.ind_Phi];
      
      // Multipoles
      Theta0_array_2D[ix-1+idx_tc_transition][ik]    = y_after_tc_solutions[ix][Constants.ind_start_theta];
      Theta1_array_2D[ix-1+idx_tc_transition][ik]    = y_after_tc_solutions[ix][Constants.ind_start_theta+1];
      Theta2_array_2D[ix-1+idx_tc_transition][ik]    = y_after_tc_solutions[ix][Constants.ind_start_theta+2];
      Theta3_array_2D[ix-1+idx_tc_transition][ik]    = y_after_tc_solutions[ix][Constants.ind_start_theta+3];
      // Derivatives
      dTheta0_dx_array_2D[ix-1+idx_tc_transition][ik]= y_deriv_after_tc_sol[ix][Constants.ind_start_theta];
      dTheta1_dx_array_2D[ix-1+idx_tc_transition][ik]= y_deriv_after_tc_sol[ix][Constants.ind_start_theta+1];
      dPhi_dx_array_2D[ix-1+idx_tc_transition][ik]   = y_deriv_after_tc_sol[ix][Constants.ind_Phi];
      dv_b_dx_array_2D[ix-1+idx_tc_transition][ik]   = y_deriv_after_tc_sol[ix][Constants.ind_vb];
      
      // Maybe implement expressions for dTheta2dx and dTheta3dx if have time in next milestone
      // for now get them from splines
      
      // Psi not dynamical so not in ODE solution.
      Psi_array_2D[ix-1+idx_tc_transition][ik]       = - Phi_array_2D[ix-1+idx_tc_transition][ik]
        - 12.0*H_0_squared/(ck_squared*a_squared)*OmegaR0*Theta2_array_2D[ix-1+idx_tc_transition][ik];
    }
    // End storing full data
    //////////////////////////////////////////////////////////////////////////////

  }
  Utils::EndTiming("integrateperturbation");

  // Check for duplicate values in the arrays, if ODE_VERBOSE_DEBUGGER in makefile
  #ifdef _FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE
  for (int ik = 0; ik < n_k; ik++)
  {
    for (int ix = 0; ix < n_x-1; ix++)
    {
      if (Theta0_array_2D[ix][ik] == Theta0_array_2D[ix+1][ik])
      {
        printf("Duplciate value in 2D array, ix: %d, Theta0(ix): %e, Theta0(ix+1): %e\n"\
          ,ix,Theta0_array_2D[ix][ik],Theta0_array_2D[ix+1][ik]);
      }
    }
  }
  #endif

  //=============================================================================
  // Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // Scalar quantities
  delta_cdm_spline.create(x_array_full,k_array,delta_cdm_array_2D,"delta_cdm_spline");
  v_cdm_spline.create(x_array_full,k_array,v_cdm_array_2D,"v_cdm_spline");
  delta_b_spline.create(x_array_full,k_array,delta_b_array_2D,"delta_b_spline");
  v_b_spline.create(x_array_full,k_array,v_b_array_2D,"v_b_spline");
  dv_b_dx_spline.create(x_array_full,k_array,v_b_array_2D,"dv_b_dx_spline");
  Phi_spline.create(x_array_full,k_array,Phi_array_2D,"Phi_spline");
  dPhi_dx_spline.create(x_array_full,k_array,dPhi_dx_array_2D,"dPhi_dx_spline");
  Psi_spline.create(x_array_full,k_array,Psi_array_2D,"Psi_spline");
  Pi_spline.create(x_array_full,k_array,Theta2_array_2D,"Pi_spline"); // Not including polarization so set to Theta2
  // Multipoles
  Theta0_spline.create(x_array_full,k_array,Theta0_array_2D,"Theta0_spline");
  Theta1_spline.create(x_array_full,k_array,Theta1_array_2D,"Theta1_spline");
  Theta2_spline.create(x_array_full,k_array,Theta2_array_2D,"Theta2_spline");
  vector_of_Theta_splines.push_back(Theta0_spline);
  vector_of_Theta_splines.push_back(Theta1_spline);
  vector_of_Theta_splines.push_back(Theta2_spline);
  // Multipole derivs
  dTheta0_dx_spline.create(x_array_full,k_array,dTheta0_dx_array_2D,"dTheta0_dx_spline");
  dTheta1_dx_spline.create(x_array_full,k_array,dTheta0_dx_array_2D,"dTheta1_dx_spline");
  vector_of_dTheta_dx_splines.push_back(Theta0_spline);
  vector_of_dTheta_dx_splines.push_back(Theta1_spline);
}

//====================================================
// Set IC at the start of the run (this is in the tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);
  
  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];

  // Set the initial conditions in the tight coupling regime
  double Psi = -2.0/3.0;

  // Scalar quantities (Phi, delta, v, ...)
  Phi       = -Psi;
  delta_cdm = -1.5*Psi;
  delta_b   = delta_cdm;
  v_cdm     = -Constants.c*k*Psi/(2.0*cosmo->Hp_of_x(x));
  v_b       = v_cdm;
  // Photon temperature perturbations (Theta_ell)
  Theta[0]  = -0.5*Psi;
  Theta[1]  = -v_cdm/3.0;
  
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
 
  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc_end_point[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc_end_point[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc_end_point[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc_end_point[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc_end_point[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc_end_point[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];
  
  // Scalar quantities (Gravitational potential, baryons and CDM)
  Phi        = Phi_tc;
  delta_b    = delta_b_tc;
  v_b        = v_b_tc;
  delta_cdm  = delta_cdm_tc;
  v_cdm      = v_cdm_tc;

  // Photon temperature perturbations (Theta_ell)
  const double ck_over_H_p_dtaudx = Constants.c*k/(cosmo->Hp_of_x(x)*rec->dtaudx_of_x(x));
  Theta[0]   = Theta_tc[0];
  Theta[1]   = Theta_tc[1];
  Theta[2]   = - 20.0/45.0*ck_over_H_p_dtaudx*Theta[1];
  for (int l = 3; l < n_ell_theta; l++)
  {
    Theta[l] = - l/(2.0*l+1.0) * ck_over_H_p_dtaudx * Theta[l-1];
  }

  return y;
}

//====================================================
// The right hand side of the perturbations ODE in the tight coupling regime
//====================================================
// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];

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
    *((OmegaCDM0*delta_cdm + OmegaB0*delta_b)/a + 4.0*OmegaR0*Theta[0]/a_squared);
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

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];

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
    *((OmegaCDM0*delta_cdm + OmegaB0*delta_b)/a + 4.0*OmegaR0*Theta[0]/a_squared);
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
  // Expression for l_max
  dThetadx[n_ell_theta-1] = ck_over_H_p*Theta[n_ell_theta-2]
  + Theta[n_ell_theta-1]*(dtaudx - Constants.c*(n_ell_theta-1+1.0)/(H_p*cosmo->eta_of_x(x)));

  return GSL_SUCCESS;
}

//====================================================
// After integrating the perturbation compute the source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  // Vector to store source function results. using 2D as in the integration solver
  Vector2D ST_array_2D(n_x,Vector(n_k));

  // Vectors to store the individual terms of the source function for debugging
  Vector2D term1(n_x,Vector(n_k));
  Vector2D term2(n_x,Vector(n_k));
  Vector2D term3(n_x,Vector(n_k));
  Vector2D term4(n_x,Vector(n_k));

  Vector2D term2_1(n_x,Vector(n_k));
  Vector2D term2_2(n_x,Vector(n_k));
  Vector2D term3_1(n_x,Vector(n_k));
  Vector2D term3_2(n_x,Vector(n_k));
  Vector2D term3_3(n_x,Vector(n_k));
  

  // Compute source functions
  for(auto ik = 0; ik < n_k; ik++){
    const double k = k_array[ik];
    const double k_Mpc = k*Constants.Mpc;
    
    for(auto ix = 0; ix < x_array_full.size(); ix++){
      const double x = x_array_full[ix];

      // Cosmological values
      const double H_p            = cosmo->Hp_of_x(x);
      const double dH_p_dx        = cosmo->dHpdx_of_x(x);
      const double ddH_p_ddx      = cosmo->ddHpddx_of_x(x);

      // Recombination values
      const double g_tilde        = rec->g_tilde_of_x(x);
      const double dg_tilde_dx    = rec->dgdx_tilde_of_x(x);
      const double ddg_tilde_ddx  = rec->ddgddx_tilde_of_x(x);
      const double tau            = rec->tau_of_x(x);

      // Perturbation values
      const double Theta0         = get_Theta(x,k,0);
      const double Psi            = get_Psi(x,k);
      const double dPsi_dx        = Psi_spline.deriv_x(x,k);
      const double dPhi_dx        = dPhi_dx_spline(x,k);
      const double Pi             = get_Pi(x,k);
      const double dPi_dx         = Pi_spline.deriv_x(x,k);
      const double ddPi_ddx       = Pi_spline.deriv_xx(x,k);
      const double v_b            = get_v_b(x,k);
      const double dv_b_dx        = dv_b_dx_spline(x,k);

      // Constants
      const double ck             = Constants.c*k;

      // Temperature source terms
      term1[ix][ik] = g_tilde*(Theta0 + Psi + 0.25*Pi);
      term2[ix][ik] = exp(-tau)*(dPsi_dx - dPhi_dx);
      term2_1[ix][ik] = dPsi_dx;
      term2_2[ix][ik] = dPhi_dx;

      term3[ix][ik] = - (dH_p_dx*g_tilde*v_b + H_p*dg_tilde_dx*v_b + H_p*g_tilde*dv_b_dx)/ck;
      term3_1[ix][ik] = dH_p_dx*g_tilde*v_b/ck;
      term3_2[ix][ik] = H_p*dg_tilde_dx*v_b/ck;
      term3_3[ix][ik] = H_p*g_tilde*dv_b_dx/ck;

      term4[ix][ik] = 3.0/(4.0*ck*ck)*(
                      Pi*g_tilde*(H_p*ddH_p_ddx + dH_p_dx*dH_p_dx)
                      + 3*H_p*dH_p_dx*(dg_tilde_dx*Pi + g_tilde*dPi_dx)
                      + H_p*H_p*(ddg_tilde_ddx*Pi + 2*dg_tilde_dx*dPi_dx + g_tilde*ddPi_ddx)
                      );
      // Temperature source
      ST_array_2D[ix][ik] = term1[ix][ik]
                            + term2[ix][ik]
                            + term3[ix][ik]
                            + term4[ix][ik];
    }
  }
  // Spline the source function
  ST_spline.create(x_array_full, k_array, ST_array_2D, "Source_Temp_x_k");

  // Debugging source terms
  std::string data_path ("../data_testing/");
  std::ofstream fp_k_values(data_path + "perturbations_k_values.txt");
  fp_k_values << "k_values per Mpc: | k_values: | horizon entry (x):\n";

  for (int ik = 0; ik < n_k; ik+=20)
  { 
    if (ik==0)
    {
      continue;
    }
    
    const double k = k_array[ik];
    const double k_Mpc = k*Constants.Mpc;
    // Find and write horizon entry
    printf("ik: %d, k: %e, k_Mpc: %.10f\n",ik,k,k_Mpc);
    double horizon_entry_x = Utils::binary_search_for_value(cosmo->eta_of_x_spline,1.0/k,
                                  std::pair(Constants.x_start,Constants.x_end_cosmo));
    fp_k_values << std::fixed << std::setprecision(5) <<  k_Mpc << " | " << 
        std::scientific << k << " | " << horizon_entry_x << "\n";

    // Configure filename and write output
    std::ostringstream stream_kvales;
    stream_kvales << std::fixed << std::setprecision(5) << k_Mpc;
    std::string filename = data_path + "testing_perturbations_k" + stream_kvales.str() + ".txt";
    std::string filename_components = data_path + "component_test_k" + stream_kvales.str() + ".txt";

    std::ofstream ST_fp(filename.c_str());
    std::cout << "Writing output to " << filename << "\n";
    std::ofstream comp_fp(filename_components.c_str());

    comp_fp << "k-value: " << k_Mpc << " " << k << "\n";
    comp_fp << std::setw(15) << "x";
    comp_fp << std::setw(15) << "dPsi";
    comp_fp << std::setw(15) << "dPhi";
    
    comp_fp << std::setw(15) << "term3_1";
    comp_fp << std::setw(15) << "term3_2";
    comp_fp << std::setw(15) << "term3_3";
    
    comp_fp << "\n";

    ST_fp << "k-value: " << k_Mpc << " " << k << "\n";
    ST_fp << std::setw(15) << "x";
    ST_fp << std::setw(15) << "ST";
    ST_fp << std::setw(15) << "ST*j_5";
    ST_fp << std::setw(15) << "ST*j_50";
    ST_fp << std::setw(15) << "ST*j_100";
    ST_fp << std::setw(15) << "arg";
    ST_fp << std::setw(15) << "term1";
    ST_fp << std::setw(15) << "term2";
    ST_fp << std::setw(15) << "term3";
    ST_fp << std::setw(15) << "term4";
    ST_fp << std::setw(15) << "j_5";
    #define ell_loop int ell = 50; ell < 251; ell+=50
    for (ell_loop)
    {
    ST_fp << std::setw(15) << "j_" + std::to_string(ell);
    }
    ST_fp << "\n";

    for (int ix = 0; ix < n_x; ix++)
    { 
      const double x = x_array_full[ix];
      if (x<0)
        {
          double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
          ST_fp << std::setw(15) << x;
          ST_fp << std::setw(15) << ST_array_2D[ix][ik];
          ST_fp << std::setw(15) << ST_array_2D[ix][ik]*Utils::j_ell(5,arg);
          ST_fp << std::setw(15) << ST_array_2D[ix][ik]*Utils::j_ell(50,arg);
          ST_fp << std::setw(15) << ST_array_2D[ix][ik]*Utils::j_ell(100,arg);
          ST_fp << std::setw(15) << arg;
          ST_fp << std::setw(15) << term1[ix][ik];
          ST_fp << std::setw(15) << term2[ix][ik];
          ST_fp << std::setw(15) << term3[ix][ik];
          ST_fp << std::setw(15) << term4[ix][ik];
          ST_fp << std::setw(15) << Utils::j_ell(5,arg);

          comp_fp << std::setw(15) << x;
          comp_fp << std::setw(15) << term2_1[ix][ik];
          comp_fp << std::setw(15) << term2_2[ix][ik];
          
          comp_fp << std::setw(15) << term3_1[ix][ik];
          comp_fp << std::setw(15) << term3_2[ix][ik];
          comp_fp << std::setw(15) << term3_3[ix][ik];

          for (ell_loop)
            {
            ST_fp << std::setw(15) << Utils::j_ell(ell,arg);
            }
          ST_fp << "\n";
          comp_fp << "\n";
        }
    }
    ST_fp.close();
    comp_fp.close();
  }
  fp_k_values.close();
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
  for (int i = 1; i < n_x; i++)
  {
    // dtaudx always negative, instead of using absolute value just force it to be positive
    tau_prime = -rec->dtaudx_of_x(x_array_full[i]);
    if (tau_prime < 10.0 || tau_prime < 10.0*k*Constants.c/cosmo->Hp_of_x(x_array_full[i]) || rec->Xe_of_x(x_array_full[i]) < 0.99)
    {
      return std::pair<double,int>(x_array_full[i-1],int(i-1)); // Return previous value and index, as condition is violated in current step
    }
  }
  // "Fail safe, to return something if the if-test fail"
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@\nDidn't find tight couple end time!\n"
    "Returning invalid x-value and index!\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  return std::pair<double,int>(999999,999999);
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
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return vector_of_Theta_splines[ell](x,k);
}

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
  fp << "k_value: " << k << "\n";
  const int npts = 20000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
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
    if (x<0)
    {
      fp << get_Source_T(x,k)  << " ";
      double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
      fp << arg << " ";
      fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
      fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
      fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
      fp << Utils::j_ell(50,  arg)           << " ";
    }
    else
    {
      fp << "0.0 0.0 0.0 0.0 0.0 0.0";
    }
    
    
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

