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
  
  Utils::StartTiming("Writing output for perturbations");
  // Experimented with different values and found these to present the different regimes
  Vector k_values = {0.3,0.1,0.013,0.007,0.001,0.0005};
  
  std::ofstream fp_k_values(data_path + "perturbations_k_values.txt");
  fp_k_values << "k_values per Mpc: | k_values: | horizon entry (x):\n";
  std::cout << "\nWriting output to " << data_path << " with following k-values per Mpc\n";
  for (auto k_value:k_values)
  {
    std::cout << k_value << std::endl;
    // Find and write horizon entry
    double horizon_entry_x = Utils::binary_search_for_value(cosmo.eta_of_x_spline,1.0/(k_value/Constants.Mpc),
                                  std::pair(Constants.x_start,Constants.x_end));
    fp_k_values << std::fixed <<  k_value << " | " << 
        std::scientific << k_value/Constants.Mpc << " | " << horizon_entry_x << "\n";

    // Configure filename and write output
    std::ostringstream stream_kvales;
    stream_kvales << std::fixed << std::setprecision(5);
    stream_kvales << k_value;
    std::string filename = "perturbations_k" + stream_kvales.str() + ".txt";
    pert.output(k_value/Constants.Mpc, data_path + filename);
  }
  fp_k_values.close();
  std::cout << std::endl;
  Utils::EndTiming("Writing output for perturbations");
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert);
  power.solve(false);
  std::string C_ell_filename = "Cells_2000k_rec_resolution.txt";
  std::string component_PS_filename = "comp_PS_2000k_rec_resolution.txt";
  std::vector<std::string> component_names{"matter"};
  power.output(data_path + C_ell_filename);
  power.output_component_power_spectrum(component_names,data_path + component_PS_filename);

  Utils::EndTiming("Everything");
}
