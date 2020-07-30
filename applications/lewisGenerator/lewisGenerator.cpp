#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include <vector>
#include <string>
#include <map>

#include <nvector/nvector_serial.h>

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

  // Compute Lewis numbers
  ComputeLewis(flame_params, flame_state_ptr);

  N_VDestroy_Serial(flame_state);

  return 0;
}
