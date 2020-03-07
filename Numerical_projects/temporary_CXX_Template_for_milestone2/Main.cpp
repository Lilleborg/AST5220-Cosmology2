#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include <string.h>

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
  // Test for input arguments, used when calculating with "toy-cosmologies"
  if (argc > 1)
  {
    // Toy cosmology used for testing in Milestone 2
    if (strcmp(argv[1],"Toy-Milestone2") == 0)
    {
      OmegaB = 0.05;
      OmegaCDM = 0.45;
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
  cosmo.output("./../data/cosmology.txt");
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

  // Output recombination quantities
  rec.output("./../data/recombination.txt");
  Utils::EndTiming("Module II");
  std::cout << "\n";
  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 * Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");
  
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
