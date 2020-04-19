#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include <string.h>

#include <sstream>
#include <iomanip>

int main(int argc, char **argv){
  std::cout << "\n";
  Utils::StartTiming("Everything");
  
  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.7;
  double OmegaB      = 0.046;
  double OmegaCDM    = 0.224;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Path for data files
  std::string data_path ("./../data/");

  // Test for input arguments, used when calculating with "toy-cosmologies"
  if (argc > 1)
  {
    // Toy cosmology used for testing
    if (strcmp(argv[1],"toy-cosmo") == 0)
    {
      OmegaB = 0.05;
      OmegaCDM = 0.45;
      data_path = "./../data_toy/";
    }
  }

  // Recombination parameters
  double Yp          = 0.0;       //Not including Helium, propper value: 0.24;

  //=========================================================================
  // Module I
  //=========================================================================
  Utils::StartTiming("Module I");
  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output(data_path + "cosmology.txt");
  Utils::EndTiming("Module I");
  std::cout << "\n";

  //=========================================================================
  // Module II
  //=========================================================================
  Utils::StartTiming("Module II");
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.info();
  rec.solve();
  rec.print_time_results();
  rec.save_time_results(data_path + "recombination_times.txt");

  // Output recombination quantities
  rec.output(data_path + "recombination.txt");
  Utils::EndTiming("Module II");
  std::cout << "\n";

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.info();
  pert.solve();
  
  Vector k_values = {0.1,0.01,0.001};

  std::ofstream fp_k_values(data_path + "perturbations_k_values.txt");
  fp_k_values << "k_values per Mpc:\n";
  for (auto k_value:k_values)
  {
    fp_k_values << k_value << "\n";

    std::ostringstream stream_kvales;
    stream_kvales << std::fixed << std::setprecision(3);
    stream_kvales << k_value;
    std::string filename = "perturbations_k" + stream_kvales.str() + ".txt";
    pert.output(k_value/Constants.Mpc, data_path + filename);
  }
  fp_k_values.close();
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert);
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
