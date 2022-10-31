#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <string>
#include <map>

#include <nvector/nvector_serial.h>

#include "utilities.h"

#include "flame_params.h"
#include "compute_lewis.h"
#include "set_initial_conditions.h"

int main(int argc, char *argv[])
{

  if(argc < 2) {
    printf("# ERROR: Incorrect command line usage.\n");
    printf("#        use %s <input parameters>\n",argv[0]);
    exit(-1);
  }


  FlameParams flame_params(argv[1]);
  N_Vector flame_state;
  flame_state = NULL;
  double *flame_state_ptr;
  flame_state_ptr = NULL;

  const int num_states = flame_params.reactor_->GetNumStates();
  int num_species = flame_params.reactor_->GetNumSpecies();

  flame_state          = N_VNew_Serial(num_states);
  flame_state_ptr      = NV_DATA_S(flame_state);

  SetConstantInlet(flame_params, flame_state_ptr);

  // ------------ BEGIN Constrained Equibrium calc  -----------//
  /*
  if(flame_params.use_equilibrium_) {
    printf("Using equilibrium composition to evaluate Lewis numbers\n");
    int k,l;
    std::vector<double> thermoCoeffs(num_species*16);
    std::vector<double> thermoCoeffs_CEQ(num_species*15);
    flame_params.mechanism_->getThermoCoeffs(&thermoCoeffs[0]);
    for(k=0; k<num_species; k++)
      for(l=0; l<15; l++) //transpose
        thermoCoeffs_CEQ[l*num_species+k] = thermoCoeffs[k*16 + l];

    int ne;
    ne = flame_params.mechanism_->getNumElements();
    std::vector<int> Ein_int(num_species*ne);
    std::vector<double> Ein_double(num_species*ne);
    flame_params.mechanism_->getSpeciesOxygenCount(&Ein_int[0]);
    flame_params.mechanism_->getSpeciesNitrogenCount(&Ein_int[2*num_species]);
    flame_params.mechanism_->getSpeciesHydrogenCount(&Ein_int[num_species]);
    if(ne>3) flame_params.mechanism_->getSpeciesCarbonCount(&Ein_int[3*num_species]);
    if(ne>4) flame_params.mechanism_->getSpeciesArgonCount(&Ein_int[4*num_species]);
    if(ne>5) flame_params.mechanism_->getSpeciesHeliumCount(&Ein_int[5*num_species]);

    for(k=0; k<num_species*ne; k++)
      Ein_double[k] = (double)Ein_int[k];

    ofstream file1;
    file1.open("CEQ-inputs");
    file1 << num_species << "\n";
    file1 << ne << "\n";
    file1 << flame_params.pressure_ << "\n";
    file1 << flame_state_ptr[num_species+1]*flame_params.ref_temperature_ << "\n";
    for(k=0; k<num_species; k++)
      file1 << 1.0/flame_params.inv_molecular_mass_[k] << " ";
    file1 << "\n";
    for(k=0; k<num_species*15; k++)
      file1 << thermoCoeffs_CEQ[k] << " ";
    file1 << "\n";
    for(k=0; k<num_species*ne; k++)
      file1 << Ein_double[k] << " ";
    file1 << "\n";
    for(k=0; k<num_species; k++)
      file1 << flame_state_ptr[k] << " ";
    file1 << "\n";
    file1.close();

    // Call CEQ
    system("/g/g90/lapointe/CEQ/eqHPfromFile.x");

    // Read equilibrium state
    std::ifstream infile("CEQ.dat");
    std::string line;
    k=0;
    double val;
    while (std::getline(infile, line))
    {
      std::istringstream iss(line);
      iss >> val;
      if(k<num_species)
        flame_state_ptr[k] = val;
      if(k==num_species)
        flame_state_ptr[k+1] = val/flame_params.ref_temperature_;

    k++;
    }
  } else {
    printf("Using inlet composition to evaluate Lewis numbers\n");
  }
  */
  // ------------ END Constrained Equibrium calc  -----------//

  // Compute Lewis numbers
  ComputeLewis(flame_params, flame_state_ptr);

  N_VDestroy_Serial(flame_state);

  return 0;
}
