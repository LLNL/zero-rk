#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string>

#include "rhs_tests.h"

int TestRateOfProgress(zerork::mechanism *mech,
                       const double temperature,
                       const double pressure,
                       const double mole_fraction[])
{
  std::string step_name;
  const int num_species   = mech->getNumSpecies();
  //const int num_reactions = mech->getNumReactions();
  const int num_steps     = mech->getNumSteps();
  int num_bad_steps=0;
  double total_concentration = pressure/(mech->getGasConstant()*temperature);
  double *concentration;
  double *net_reaction_rate, *creation_rate, *destruction_rate, *step_rate;

  concentration     = new double[num_species];
  net_reaction_rate = new double[num_species];
  creation_rate     = new double[num_species];
  destruction_rate  = new double[num_species];
  step_rate         = new double[num_steps];

  for(int j=0; j<num_species; ++j) {
    concentration[j] = total_concentration*mole_fraction[j];
  }
  mech->getReactionRates(temperature,
                         concentration,
                         net_reaction_rate,
                         creation_rate,
                         destruction_rate,
                         step_rate);   

  for(int j=0; j<num_steps; ++j) {

    if(isnan(step_rate[j])) {
      step_name.clear();
      mech->getReactionNameDirOfStep(j,&step_name);
      printf("# NaN rate of progress (%g) for step %d: %s\n",
             step_rate[j],j,step_name.c_str());
      fflush(stdout);
      ++num_bad_steps;
    }
    else if(isinf(step_rate[j])) {
      step_name.clear();
      mech->getReactionNameDirOfStep(j,&step_name);
      printf("# Inf rate of progress (%g) for step %d: %s\n",
             step_rate[j],j,step_name.c_str());
      fflush(stdout);
      ++num_bad_steps;
    }
    // else if(isnormal(step_rate[j]) == 0) {
    //   step_name.clear();
    //   mech->getReactionNameDirOfStep(j,&step_name);
    //   printf("# Rate of progress fails isnormal(%g) for step %d: %s\n",
    //          step_rate[j],j,step_name.c_str());
    //   fflush(stdout);
    //   ++num_bad_steps;
    // }

  }  
 
  delete [] concentration;
  delete [] net_reaction_rate;
  delete [] creation_rate;
  delete [] destruction_rate;
  delete [] step_rate;

  return num_bad_steps;
}
