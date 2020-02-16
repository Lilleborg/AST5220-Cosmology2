#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM,
    double Neff, 
    double TCMB) :
  // This initializes the class variables h,OmegaB,OmegaCDM... as the value they have where sent into the class
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
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
  OmegaK = - K * Constants.c * Constants.c / H0 / H0;

  //OmegaLambda
  OmegaLambda = 1 - (OmegaR + OmegaNu + OmegaB + OmegaCDM + OmegaK);
  
  // Check sum of the Omegas today
  Omega_summed = OmegaB + OmegaCDM + OmegaLambda + OmegaR + OmegaK + OmegaNu;
  std::cout << "Sum of the present day capital omegas: " << Omega_summed << std::endl;
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
// Helper methods (private)
//====================================================

  // Expression used in all get_Omegas methods, returns (H0/H)^2
double BackgroundCosmology::H0_over_H_squared(double x) const{
  double H_temp = H_of_x(x);
  return H0*H0/H_temp/H_temp;
}

  // Returns a vector with components exp(3*x) and exp(4*x)
Vector BackgroundCosmology::exp_of_3x_and_4x(double x) const{
  double a = exp(x);
  double exp3x = a*a*a;
  double exp4x = exp3x*a;
  Vector res{exp3x,exp4x};
  return res;
}

//====================================================
// Get methods
//====================================================

  // Returns the Hubble parameter as function of x, using Friedmann 1
double BackgroundCosmology::H_of_x(double x) const{
  Vector exponentials = exp_of_3x_and_4x(x);
  double res = H0 * sqrt(
    (OmegaB+OmegaCDM) / exponentials[0]
    + OmegaR / exponentials[1]
    + OmegaLambda);

  return res;
}

  // Returns the derivative of Hubble parameter wrt x
double BackgroundCosmology::dHdx_of_x(double x) const{
  Vector exponentials = exp_of_3x_and_4x(x);
  double res = H0*H0/(2*H_of_x(x)) * (-3*(OmegaB+OmegaCDM)/exponentials[0] - 4*OmegaR/exponentials[1]);

  return res;
}

  // Returns the double derivative of Hubble parameter wrt x
double BackgroundCosmology::ddHddx_of_x(double x) const{
  double H = H_of_x(x);
  Vector exponentials = exp_of_3x_and_4x(x);
  double res = H0*H0/2 * (1/H * (9*(OmegaB+OmegaCDM)/exponentials[0] + 16*OmegaR/exponentials[1]) - dHdx_of_x(x)/H/H * (-3*(OmegaB+OmegaCDM)/exponentials[0] - 4*OmegaR/exponentials[1]));
  return res;
}
  // Returns Hubble prime as function of x and H_of_x
double BackgroundCosmology::Hp_of_x(double x) const{
  double res = exp(x) * H_of_x(x);

  return res;
}


  // Returns the derivative of Hubble prime wrt x
double BackgroundCosmology::dHpdx_of_x(double x) const{
  double a = exp(x);
  double res = a*dHdx_of_x(x) + H_of_x(x)*a;

  return res;
}

  // Returns the double derivative of Hubble prime wrt x
double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double a = exp(x);
  double dHdx = dHdx_of_x(x);
  double H = H_of_x(x);
  double ddHddx = ddHddx_of_x(x);

  double res = a*(2*dHdx + H + ddHddx);

  return res;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  Vector exponentials = exp_of_3x_and_4x(x);
  double Omega = H0_over_H_squared(x) * OmegaB / exponentials[0];

  return Omega;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  Vector exponentials = exp_of_3x_and_4x(x);
  double Omega = H0_over_H_squared(x) * OmegaR / exponentials[1];

  return Omega;
}

  // Calculating without neutrinos, so always returning zero
double BackgroundCosmology::get_OmegaNu(double x) const{ 
  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  Vector exponentials = exp_of_3x_and_4x(x);
  double Omega = H0_over_H_squared(x) * OmegaCDM / exponentials[0];

  return Omega;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  double Omega = H0_over_H_squared(x) * OmegaLambda;

  return Omega;
}

  // Calculating without curvature, so always returning zero
double BackgroundCosmology::get_OmegaK(double x) const{
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
  const double x_max = 2.0;
  const int    n_pts = npts;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";

    fp << H_of_x(x)          << " ";
    fp << dHdx_of_x(x)       << " ";

    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";

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

