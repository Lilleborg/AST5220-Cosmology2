#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaLambda,
    double Neff, 
    double TCMB) :
  // This initializes the class variables h,OmegaB,OmegaCDM... as the value they have where sent into the class
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaLambda(OmegaLambda),
  Neff(Neff), 
  TCMB(TCMB)
{
  // H0:
  H0 = Constants.H0_over_h * h;
  
  //OmegaR:
  OmegaR = M_PI*M_PI / 15 * pow(Constants.k_b * TCMB,4) / pow(Constants.hbar,3) / pow(Constants.c,5)
    * 8 * M_PI * Constants.G / 3 / H0 / H0;
  
  //OmegaNu:
  OmegaNu = 0;

  //OmegaK:
  const double K = 0;
  OmegaK = K / H0 / H0;
  
  // Check sum of the Omegas today
  Omega_summed = OmegaB + OmegaCDM + OmegaLambda + OmegaR + OmegaK + OmegaNu;
  if (Omega_summed != 1.0) // stricked test for unity
    {
    std::cout << "Sum of the present day capital omegas, (should be unity): " << Omega_summed << std::endl;
    }
}

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  Vector x_array = Utils::linspace(x_start,x_end,npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c / Hp_of_x(x);
    return GSL_SUCCESS;
  };

  // Solving the ODE:
  Vector eta_ini{0.0};
  ODESolver ode;
  ode.solve(detadx,x_array,eta_ini);
  auto eta_all_data = ode.get_data();
  // Filling data into array:
  Vector eta_array(x_array.size());
  for (int i = 0; i < eta_all_data.size(); i++){
    eta_array.at(i) = eta_all_data[i][0];
  }

  // Creating spline:
  eta_of_x_spline.create(x_array,eta_array);
  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  // Returns the Hubble parameter as function of x, using Friedmann 1
  double H_temp = H0 * sqrt(
    (OmegaB+OmegaCDM) * exp(-3 * x)
    + OmegaR * exp(-4 * x)
    + OmegaLambda);

  return H_temp;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  // Returns Hubble prime as function of x and H_of_x
  double Hprime_temp = exp(x) * H_of_x(x);

  return Hprime_temp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::H0_over_H_squared(double x) const{
  double H_temp = H_of_x(x);
  return H0*H0/H_temp/H_temp;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  double Omega = H0_over_H_squared(x) * OmegaB * exp(-3 * x);

  return Omega;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  double Omega = H0_over_H_squared(x) * OmegaR * exp(-4 * x);

  return Omega;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  // Calculating without neutrinos, so always returning zero
  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  double Omega = H0_over_H_squared(x) * OmegaCDM * exp(-3 * x);

  return Omega;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  double Omega = H0_over_H_squared(x) * OmegaLambda;

  return Omega;
}

double BackgroundCosmology::get_OmegaK(double x) const{
  // Calculating without curvature, so always returning zero
  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB() const{ 
  return TCMB; 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -15;
  const double x_max = 0.0;
  const int    n_pts = npts;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

