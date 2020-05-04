#include "Utils.h"

struct ConstantsAndUnits Constants;

namespace Utils{

  static std::map<std::string, std::chrono::time_point<std::chrono::steady_clock>> timings;
  
  void StartTiming(const char *name){
    std::string sname(name);
    StartTiming(sname);
  }
  void EndTiming(const char *name){
    std::string sname(name);
    EndTiming(sname);
  }

  void StartTiming(std::string &name){
    std::chrono::time_point<std::chrono::steady_clock> start_time = std::chrono::steady_clock::now();
    timings[name] = start_time;
  }
  void EndTiming(std::string &name){
    std::chrono::time_point<std::chrono::steady_clock> end_time = std::chrono::steady_clock::now();
   
    auto it = timings.find(name);
    if (it == timings.end()){
      std::cout << "Can't print elapsed time. Start time for time point [" << name << "] was not found\n"; 
    } else {
      auto start_time = it->second;
      std::cout << "Elapsed time for [" << name << "]: " << timeInSeconds(start_time, end_time) << " sec\n";
      timings.erase(it);
    }
  }

  // For timing
  std::chrono::time_point<std::chrono::steady_clock> getTime(){
    return std::chrono::steady_clock::now();
  }

  double timeInSeconds(
      std::chrono::time_point<std::chrono::steady_clock> & time_start, 
      std::chrono::time_point<std::chrono::steady_clock> & time_end){
    return std::chrono::duration_cast<std::chrono::duration<double>> (time_end - time_start).count();
  }

  // Find a value in a spline
  double binary_search_for_value(
      Spline const &y, 
      double y_value, 
      std::pair<double,double> xrange, 
      double epsilon){
    const int nsearch   = 20;
    const int nmax_iter = 1000;

    // The search range. If not provided use the full range from the spline
    if(xrange.first == 0.0 && xrange.second == 0.0)
      xrange = y.get_xrange();

    // Arrange search ranghe so that x_low < x_high)
    double x_low   = std::min(xrange.first, xrange.second);
    double x_high  = std::max(xrange.first, xrange.second);

    // Accuracy we want
    double delta_x = (x_high - x_low) * epsilon; 

    // We measure function values relative to y_value
    double y_high = y(x_high) - y_value;
    double y_low  = y(x_low)  - y_value;

    // If the endpoints are both larger than zero do 
    // a coarse grid search to see if we can find a good value
    // If not then throw an error
    if(y_low * y_high >= 0.0){
      // Make a search grid
      auto x_array = linspace(x_low, x_high, nsearch);

      // Are the points positive or negative?
      int sign = y_low > 0.0 ? 1 : -1;

      for(int i = 0; i < nsearch; i++){
        double ycur = y(x_array[i]) - y_value;
        if(sign * ycur < 0.0){
          y_low = ycur;
          x_low = x_array[i];
          break;
        }
        if(i == nsearch-1) throw "Error binary_search_for_value. Could not find a good interval to start search\n";
      }
    }

    // Make the function we are going to evaluate adjusted so that y_low < y_high
    int sign = y_high > y_low ? 1 : -1;
    auto function = [&](double x){
      return (y(x) - y_value) * sign;
    };

    // Do the binary search
    int count = 0;
    while(x_high - x_low > delta_x){
      double x_mid = (x_high+x_low)/2.0;
      double y_mid = function(x_mid);
      if(y_mid > 0.0)
        x_high = x_mid;
      else
        x_low  = x_mid;
      count++;
      if(count > nmax_iter) throw "Error binary_search_for_value. Value not found after nmax iterations\n";
    }

    return x_low;
  }

  // Compute j_ell(x) for all ell = 0, ..., lmax using a recursion function with input of exact values 
  std::vector<double> j_ell_array(int lmax, const double x){
    std::vector<double> res(lmax+1, 0.0);

    // We need some input
    const double epsilon = 1e-100;
    const double jellmax       = j_ell(lmax  ,x);
    const double jellmaxplus1  = j_ell(lmax+1,x);
    const double jellzero      = j_ell(0     ,x);

    // Recursion relation for j_(n+1) / jn
    double h = fabs(jellmax) < epsilon ? 0.0 : jellmaxplus1 / jellmax;
    for(int k = lmax; k >= 1; k--){
      h = x/(2*k+1 - x * h);
      res[k] = h;
    }

    // Transform into j_n
    res[0] = jellzero;
    for(int ell = 1; ell <= lmax; ell++){
      res[ell] *= res[ell-1];
    }

    return res;
  }

  // Useful function for generating a linspace
  std::vector<double> linspace(double xmin, double xmax, int num){
    std::vector<double> res(num);
    double delta_x = (xmax-xmin)/double(num-1);
    for(int i = 0; i < num; i++){
      res[i] = xmin + delta_x * i;
    }
    res[num-1] = xmax; // Just to make sure its exactly xmax!
    return res;
  }

  // Function to generate a logspace with default base e
  std::vector<double> logspace(double min_exponent, double max_exponent, int num, double base){
    std::vector<double> res(num);
    std::vector<double> exponents = linspace(min_exponent,max_exponent,num);
    for (int i = 0; i < num; i++)
    {
      res[i] = pow(base,exponents[i]); // using pow here.. not ment to be called often
    }
    return res;
  }

  // Function to get the spherical Bessel function j_n(x)
  double j_ell(const int ell, const double x){
    if(x==0.0) return ell == 0 ? 1.0 : 0.0;

    // In this regime the function is <~ 1e-6 times the maximum value so safe to put it to zero
    // to avoid issues with the library functions failing to compute it
    if(ell >= 10 && x < (1.0 - 2.6/sqrt(ell)) * ell) return 0.0;

  #if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || (__cplusplus >= 201703L))
      // If you have a c++17 compiler you can use this
    
      // The library fails for the largest arguments so simply put to zero
      // to avoid any issues with this (these large values not very relevant for us anyway)
      if(x > 14000.0) return 0.0;
    
      return std::sph_bessel(ell, x);
  #else
      // Otherwise lets use GSL 
   
      // For 'small' ell GSL fails for the largest arguments so simply put to zero
      // to avoid issues with this (these large values not very relevant for us anyway)
      if(ell < 500 && x > 9000) return 0.0; 

      return gsl_sf_bessel_jl(ell, x);
  #endif
    }


  
  double J_n(const int n, const double x){
#ifdef _COMPLEX_BESSEL
    return sp_bessel::besselJ(n, x).real();
#else
    return gsl_sf_bessel_Jn(n, x);
#endif
  }

  // if(n > 100 && x < 0.2 * n) return 0.0;    
  //  return gsl_sf_bessel_Jn(n, x); 

  // // Function to get the spherical Bessel function j_n(x)
  //   double j_ell(const int n, const double x){
  //     if(x==0.0) return n == 0 ? 1.0 : 0.0;

  // #if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || (__cplusplus >= 201703L))
  //     // If you have a c++17 compiler you can use this
  //     return std::sph_bessel(n, x);
  // #else
  //     // Otherwise lets use GSL with approximation for x << n to prevent 
  //     // underflow issues in that implementation for large n and small x
  //     // std::cout << "Using GSL \n";
  //     if(n > 100 && x < 0.2 * n) return 0.0;
  //     return gsl_sf_bessel_jl(n, x);
  // #endif
  //   }
  //   double J_n(const int n, const double x){
  // #if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || (__cplusplus >= 201703L))
  //     // If you have a c++17 compiler you can use this
  //     return std::cyl_bessel(n, x);
  // #else
  //     if(n > 100 && x < 0.2 * n) return 0.0;
  //     return gsl_sf_bessel_Jn(n, x);
  // #endif
  //   }

  // 2-point stencil with zero derivative at the end-points
  std::vector<double> derivative(std::vector<double> &x, std::vector<double> &f){
    int n = f.size();
    std::vector<double> derivative(n);
    for(int i = 1; i < n-1; i++){
      double dx   = x[i+1] - x[i];
      derivative[i] = (f[i+1] - f[i])   / dx;
    }
    derivative[n-1] = derivative[n-2];
    return derivative;
  }
}

// Be able to take common mathematical functions on a vector
#define FUNS(FUN) \
  struct op_##FUN { double operator() (double d) const { return FUN(d); } }; \
  std::vector<double> FUN(const std::vector<double>& x) { \
    const int n = x.size(); \
    std::vector<double> y(n); \
    std::transform(x.begin(), x.end(), y.begin(), op_##FUN()); \
    return y; \
  } 
FUNS(exp); FUNS(log); FUNS(cos); FUNS(sin); FUNS(tan); FUNS(fabs); FUNS(atan)
#undef FUNS
