#include <math.h>
#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.

#ifdef SUNDIALS2
#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR
#elif SUNDIALS3
#include <cvode/cvode_spils.h>
//#include <sunlinsol/sunlinsol_spgmr.h>
#elif SUNDIALS4
//#include <sunlinsol/sunlinsol_spgmr.h>
//#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif


#include "idtControlParams.h"
#include "idtSolvers.h" // for CVReactorState * getState(...)
#include "jacobian_stats.h"


// Comparison function in the format to be used by qsort() from stdlib.h.
// The comparison function returns an integer less than, equal to,  or
// greater  than  zero  if  the first argument is considered to be respec‐
// tively less than, equal to, or greater than the second. Specifically, the
// size of the argument is measured as the absolute value of the rate of
// progress derivative.
//
// This function, when called within, qsort will sort the array in ASCENDING
// order.
int compareROPDerivative(const void *aptr, const void *bptr)
{
  ROPDerivative a = *((ROPDerivative *)aptr);
  ROPDerivative b = *((ROPDerivative *)bptr);

  if(fabs(a.value) < fabs(b.value)) {
    return -1;
  } else if(fabs(a.value) > fabs(b.value)) {
    return 1;
  }
  return 0;
}
// Comparison function in the format to be used by qsort() from stdlib.h.
// The comparison function returns an integer greater than, equal to,  or
// less  than  zero  if  the first argument is considered to be respec‐
// tively less than, equal to, or greater than the second. Specifically, the
// size of the argument is measured as the absolute value of the rate of
// progress derivative.
//
// This function, when called within, qsort will sort the array in DESCENDING
// order.
int compareROPDerivativeDescending(const void *aptr, const void *bptr)
{return -compareROPDerivative(aptr,bptr);}

//
int initializeROPDerivativeVector(idtControlParams *idtCtrl,
                                  std::vector<ROPDerivative> *vec)
{
  int reaction_id;
  ROPDerivative rop_deriv;
  vec->clear();

  for(int j=0; j<idtCtrl->mech->getNumSteps(); ++j) {
    reaction_id = idtCtrl->mech->getRxnIdxOfStep(j);
    rop_deriv.value = 0.0;
    rop_deriv.step_id = j;
    rop_deriv.elementary_step_order = idtCtrl->mech->getOrderOfStep(j);
    rop_deriv.is_pressure_dependent = false;
    if(idtCtrl->mech->isFalloffReaction(reaction_id) == 1 ||
       idtCtrl->mech->isThirdBodyReaction(reaction_id) == 1) {
      rop_deriv.is_pressure_dependent = true;
    }
    for(int k=0; k<idtCtrl->mech->getOrderOfStep(j); ++k) {
      rop_deriv.species_id = idtCtrl->mech->getSpecIdxOfStepReactant(j,k);
      vec->push_back(rop_deriv);
    }
  }

  return static_cast<int>(vec->size());
}

double updateROPDerivativeVector(const double rate_of_progress[],
                                 const double concentration[],
                                 std::vector<ROPDerivative> *vec)
{
  double max_value=0.0;
  for(unsigned int j=0; j<vec->size(); ++j) {
    (vec->at(j)).value = rate_of_progress[(vec->at(j)).step_id]/
      concentration[(vec->at(j)).species_id];
    if(fabs(vec->at(j).value) > max_value) {
      max_value = fabs(vec->at(j).value);
    }
  }
  return max_value;
}

void sortROPDerivativeVector(std::vector<ROPDerivative> *vec)
{
  qsort((void *)vec->data(),
        vec->size(),
        sizeof(ROPDerivative),
        compareROPDerivativeDescending);
}

void getROPDerivativeVector(idtControlParams *idtCtrl,
                            std::vector<ROPDerivative> *vec,
                            double *min_concentration)
{
  double total_concentration;
  std::vector<double> net_rate;
  std::vector<double> creation_rate;
  std::vector<double> destruction_rate;
  std::vector<double> step_rate_of_progress;

  CVReactorState thermo_state;


  // clears the vector and assigns the species and step indexes of each
  // rate of progress derivative (with respect to concentration)
  initializeROPDerivativeVector(idtCtrl,vec);

  // use the assign() vector member function to make sure there is contiguous
  // memory allocated for each vector so that pointers to their data can be
  // passed to zerork functions as double[].
  net_rate.assign(idtCtrl->mech->getNumSpecies(),0.0);
  creation_rate.assign(idtCtrl->mech->getNumSpecies(),0.0);
  destruction_rate.assign(idtCtrl->mech->getNumSpecies(),0.0);
  step_rate_of_progress.assign(idtCtrl->mech->getNumSteps(),0.0);

  // get the thermodynamic state ( concentration and temperature)
  getState(idtCtrl,&thermo_state);

  // If the user
  if(*min_concentration <= 0.0) {
    total_concentration = thermo_state.pressure/
      (idtCtrl->mech->getGasConstant()*thermo_state.temperature);
    (*min_concentration) = total_concentration*idtCtrl->cvodeCtrl.absTol*
       idtCtrl->cvodeCtrl.relTol;
  }
  // set all the concentrations c[i] <= max(c[i],min_concentration)
  // NOTE: the threshold in the jacobian matrix construction is applied to the
  // mass fractions (idtCtrl->odeUserParams.minMassFrac) so the terms may not
  // match exactly.
  for(int j=0; j<idtCtrl->mech->getNumSpecies(); ++j) {
    if(thermo_state.concentration[j] < (*min_concentration)) {
      thermo_state.concentration[j] = (*min_concentration);
    }
  }

  // get the rate of progress for each step
  idtCtrl->mech->getReactionRatesLimiter(thermo_state.temperature,
                                  thermo_state.concentration.data(),
				  idtCtrl->odeUserParams.step_limiter,
                                  net_rate.data(),
                                  creation_rate.data(),
                                  destruction_rate.data(),
                                  step_rate_of_progress.data());
  // compute the ROPDerivative values
  updateROPDerivativeVector(step_rate_of_progress.data(),
                            thermo_state.concentration.data(),
                            vec);

  // sort the rate of progress derivative values by magnitude, largest first
  sortROPDerivativeVector(vec);
}

int getROPDerivative(idtControlParams *idtCtrl,
                     double *min_concentration,
                     double rate_of_progress[],
                     int species_id[],
                     int step_id[])
{
  CVReactorState thermo_state;
  std::vector<ROPDerivative> vec;

  // if the output arrays have not been allocate then just return the
  // size need for each output array
  if(rate_of_progress == NULL || species_id == NULL || step_id == NULL) {
    return idtCtrl->mech->getTotalReactants();
  }
  getROPDerivativeVector(idtCtrl,&vec,min_concentration);
  if(static_cast<int>(vec.size()) != idtCtrl->mech->getTotalReactants()) {
    printf("ERROR: In getROPDerivative(...),\n");
    printf("       getROPDerivativeVector(...) vector size %lu does not\n",
           vec.size());
    printf("       the total number of step reactants %d in the mechanism.\n",
           idtCtrl->mech->getTotalReactants());
    exit(-1);
  }
  for(unsigned int j=0; j<vec.size(); ++j) {
    species_id[j]       = vec[j].species_id;
    step_id[j]          = vec[j].step_id;
    rate_of_progress[j] = vec[j].value;
  }
  return idtCtrl->mech->getTotalReactants();
}

void writeROPDerivativeReport(idtControlParams *idtCtrl,
                              const double t_current,
                              const double min_value,
                              FILE *fptr)
{
  double time_step;
  double min_concentration=-1.0e300;
  std::vector<ROPDerivative> vec;
  long int num_jac;
#if defined SUNDIALS2 || defined SUNDIALS3
  CVSpilsGetNumPrecEvals(idtCtrl->cvodeCtrl.cvodeMemPtr,&num_jac);
#else
  CVodeGetNumPrecEvals(idtCtrl->cvodeCtrl.cvodeMemPtr,&num_jac);
#endif
  getROPDerivativeVector(idtCtrl,
                         &vec,
                         &min_concentration); // min_concentration will be
                                              // recomputed
  fprintf(fptr,
          "# Rate of Progress Derivative Report at t = %.18g [s]\n",
          t_current);

  CVodeGetLastStep(idtCtrl->cvodeCtrl.cvodeMemPtr,&time_step);
  fprintf(fptr,
          "# Last time step -    CVodeGetLastStep()    = %.18g [s]\n",
          time_step);
  CVodeGetCurrentStep(idtCtrl->cvodeCtrl.cvodeMemPtr,&time_step);
  fprintf(fptr,
          "# Current time step - CVodeGetCurrentStep() = %.18g [s]\n",
          time_step);
  fprintf(fptr,
          "# Number of Jacobians - CVSpilsGetNumPrecEvals() = %d\n",
          num_jac);
  fprintf(fptr,
          "# The concentration derivative of each step's rate of progress is sorted in\n# descending order by magnitude.  To facilitate calculation, the derivatives\n# are calculated for a strictly positive species concentration composition.\n# Specifically, the concentration C[i] of species 'i' is set as follows:\n#   C[i] = max(C_min,C[i]) with C_min = %.18g [kmol/m^3].\n",
          min_concentration);
  fprintf(fptr,
          "#\n# Note that ONLY rate of progress derivatives with a magnitude\n");
  fprintf(fptr,"# greater than %12.3e [Hz] are reported.\n",min_value);
  fprintf(fptr,
          "#\n#       [Hz]   species  step\n");
  fprintf(fptr,
          "#    d(ROP)/dC, index, index, description\n");
  for(unsigned int j=0; j<vec.size(); ++j) {
    if(fabs(vec[j].value) > min_value) {
      int reaction_id = idtCtrl->mech->getRxnIdxOfStep(vec[j].step_id);
      char reaction_dir[]="(fwd)";
      if(idtCtrl->mech->getStepIdxOfRxn(reaction_id,-1) == vec[j].step_id) {
        strcpy(reaction_dir,"(rev)");
      }
      fprintf(fptr,
              "%14.7e %6d %6d  d/dC[%s] { %s %s }\n",
              vec[j].value,
              vec[j].species_id,
              vec[j].step_id,
              idtCtrl->mech->getSpeciesName(vec[j].species_id),
              idtCtrl->mech->getReactionName(reaction_id),
              reaction_dir);
    }
  }
}
// Calculates the rate-of-progress temperature derivative using the
// finite difference approximation.  The derivatives are stored in the
// ROPDerivative vector and are NOT sorted by magnitude.  If min_concentration
// is set to NULL, then there is no strictly positive correction to the
// concentration composition.
void getROPTempDerivativeVector(idtControlParams *idtCtrl,
                               std::vector<ROPDerivative> *vec,
                               double *min_concentration)
{
  double multiplier;
  double total_concentration;
  double temperature2, delta_temperature;
  std::vector<double> net_rate;
  std::vector<double> creation_rate;
  std::vector<double> destruction_rate;
  std::vector<double> step_rate_of_progress1,step_rate_of_progress2;

  CVReactorState thermo_state;
  ROPDerivative initial;

  // clears the vector and assigns the species and step indexes of each
  // rate of progress derivative (with respect to concentration)
  vec->clear();

  for(int j=0; j<idtCtrl->mech->getNumSteps(); ++j) {
    initial.species_id = temperature_id;
    initial.step_id = j;
    initial.value = 0.0;
    vec->push_back(initial);
  }

  // use the assign() vector member function to make sure there is contiguous
  // memory allocated for each vector so that pointers to their data can be
  // passed to zerork functions as double[].
  net_rate.assign(idtCtrl->mech->getNumSpecies(),0.0);
  creation_rate.assign(idtCtrl->mech->getNumSpecies(),0.0);
  destruction_rate.assign(idtCtrl->mech->getNumSpecies(),0.0);
  step_rate_of_progress1.assign(idtCtrl->mech->getNumSteps(),0.0);
  step_rate_of_progress2.assign(idtCtrl->mech->getNumSteps(),0.0);

  // get the thermodynamic state ( concentration and temperature)
  getState(idtCtrl,&thermo_state);
  total_concentration = thermo_state.pressure/
    (idtCtrl->mech->getGasConstant()*thermo_state.temperature);

  if(min_concentration != NULL) {
    // only apply the minimum concentration if not set to NULL
    if(*min_concentration <= 0.0) {
      (*min_concentration) = total_concentration*idtCtrl->cvodeCtrl.absTol*
         idtCtrl->cvodeCtrl.relTol;
    }
    // set all the concentrations c[i] <= max(c[i],min_concentration)
    // NOTE: the threshold in the jacobian matrix construction is applied to the
    // mass fractions (idtCtrl->odeUserParams.minMassFrac) so the terms may not
    // match exactly.
    for(int j=0; j<idtCtrl->mech->getNumSpecies(); ++j) {
      if(thermo_state.concentration[j] < (*min_concentration)) {
        thermo_state.concentration[j] = (*min_concentration);
      }
    }
  }

  // get the rate of progress for each step for the original temperature
  idtCtrl->mech->getReactionRatesLimiter(thermo_state.temperature,
                                  thermo_state.concentration.data(),
				  idtCtrl->odeUserParams.step_limiter,
                                  net_rate.data(),
                                  creation_rate.data(),
                                  destruction_rate.data(),
                                  step_rate_of_progress1.data());

  delta_temperature = thermo_state.temperature*
    idtCtrl->odeUserParams.sqrtUnitRnd;
  temperature2 = thermo_state.temperature+delta_temperature;
  delta_temperature = temperature2-thermo_state.temperature;
  // get the rate of progress for each step for the original temperature
  idtCtrl->mech->getReactionRatesLimiter(temperature2,
                                  thermo_state.concentration.data(),
				  idtCtrl->odeUserParams.step_limiter,
                                  net_rate.data(),
                                  creation_rate.data(),
                                  destruction_rate.data(),
                                  step_rate_of_progress2.data());

  // compute the ROP Temperature Derivative values
  multiplier = idtCtrl->refTemp/(delta_temperature*total_concentration);
  for(int j=0; j<idtCtrl->mech->getNumSteps(); ++j) {
    (*vec)[j].value = multiplier*(step_rate_of_progress2[j] -
                                  step_rate_of_progress1[j]);
  }
}

void writeROPTempDerivativeReport(idtControlParams *idtCtrl,
                                  const double t_current,
                                  const double min_value,
                                  FILE *fptr)
{
  double min_concentration=-1.0e300;
  double total_concentration;
  double norm_factor;
  CVReactorState thermo_state;
  std::vector<ROPDerivative> vec,vec_nocorrection;

  getState(idtCtrl,&thermo_state);
  total_concentration = thermo_state.pressure/
      (idtCtrl->mech->getGasConstant()*thermo_state.temperature);

  getROPTempDerivativeVector(idtCtrl,
                             &vec,
                             &min_concentration); // min_concentration will be
                                                  // recomputed
  getROPTempDerivativeVector(idtCtrl,
                             &vec_nocorrection,
                             NULL); // do not apply the strictly positive
                                    // correction to the concentrations
  sortROPDerivativeVector(&vec_nocorrection);
  fprintf(fptr,
          "# Rate of Progress Temperature Derivative Report at t = %.18g [s]\n",
          t_current);
  fprintf(fptr,
          "# The temperature derivative of each step's rate of progress is sorted in\n# descending order by magnitude for the uncorrected concentration.  To\n# facilitate calculation, the derivatives are calculated using a finite\n# difference approximation. The temperature derivatives using a strictly for\n# a strictly positive species concentration composition are also given.\n# Specifically, the concentration C[i] of species 'i' is set as follows:\n#   C[i] = max(C_min,C[i]) with C_min = %.18g [kmol/m^3].\n",
          min_concentration);
  norm_factor = idtCtrl->refTemp/total_concentration;
  fprintf(fptr,"#   A_ref = T_ref/C_tot = %.18g [K*m^3/kmol]\n",
          norm_factor);

  fprintf(fptr,"#     with T_ref = %.18g [K]\n",idtCtrl->refTemp);
  fprintf(fptr,"#     and  C_tot = %.18g [kmol/m^3]\n",total_concentration);
  fprintf(fptr,
          "#\n# Note that ONLY rate of progress derivatives with a magnitude\n");
  fprintf(fptr,"# greater than %12.3e [Hz] are reported.\n",min_value);

  fprintf(fptr,
          "#\n#       [Hz]             [Hz]\n");
  fprintf(fptr,
          "#     uncorrected      C[i]>=C_min    step\n");
  fprintf(fptr,
          "#   A_ref*d(ROP)/dT, A_ref*d(ROP)/dT, index, description\n");

  for(unsigned int j=0; j<vec_nocorrection.size(); ++j) {

    if(fabs(vec_nocorrection[j].value) > min_value) {
      char reaction_dir[]="(fwd)";
      int reaction_id =
        idtCtrl->mech->getRxnIdxOfStep(vec_nocorrection[j].step_id);

      if(idtCtrl->mech->getStepIdxOfRxn(reaction_id,-1) ==
        vec_nocorrection[j].step_id) {
          strcpy(reaction_dir,"(rev)");
      }

      fprintf(fptr,
              "     %14.7e   %14.7e %6d  A_ref*d/dT { %s %s }\n",
              vec_nocorrection[j].value,
              vec[vec_nocorrection[j].step_id].value,
              vec_nocorrection[j].step_id,
              idtCtrl->mech->getReactionName(reaction_id),
              reaction_dir);
    }
  }
}
void updateROPDerivativeDistribution(idtControlParams *idtCtrl,
                                     const double t_current,
                                     zerork::utilities::Distribution *order1_cum,
                                     zerork::utilities::Distribution *order2_cum,
                                     zerork::utilities::Distribution *order3g_cum,
                                     zerork::utilities::Distribution *thirdfall_cum,
                                     zerork::utilities::Distribution *temp_cum,
                                     FILE *fptr)
{
  long int num_jac;
  double dt_last;
  double min_concentration = -1.0e300;
  std::vector<ROPDerivative> d_concentration, d_temperature;
  zerork::utilities::Distribution order1_dist=*order1_cum;
  zerork::utilities::Distribution order2_dist=*order2_cum;
  zerork::utilities::Distribution order3g_dist=*order3g_cum;
  zerork::utilities::Distribution thirdfall_dist=*thirdfall_cum;
  zerork::utilities::Distribution temp_dist=*temp_cum;
  CVReactorState thermo_state;

  order1_dist.zeroBinWeights();
  order2_dist.zeroBinWeights();
  order3g_dist.zeroBinWeights();
  thirdfall_dist.zeroBinWeights();
  temp_dist.zeroBinWeights();

  // get the thermodynamic state ( concentration and temperature)
  getState(idtCtrl,&thermo_state);

  getROPDerivativeVector(idtCtrl,
                         &d_concentration,
                         &min_concentration);
  getROPTempDerivativeVector(idtCtrl,
                             &d_temperature,
                             &min_concentration);
  // add ROP derivatives
  for(unsigned int j=0; j<d_concentration.size(); ++j) {
    int step_id = d_concentration[j].step_id;
    int reaction_id = idtCtrl->mech->getRxnIdxOfStep(step_id);
    if(idtCtrl->mech->isThirdBodyReaction(reaction_id) == 1 ||
       idtCtrl->mech->isFalloffReaction(reaction_id) == 1) {

      thirdfall_dist.addValue(d_concentration[j].value,1.0);

    } else if(idtCtrl->mech->getOrderOfStep(step_id) == 1) {

      order1_dist.addValue(d_concentration[j].value,1.0);

    } else if(idtCtrl->mech->getOrderOfStep(step_id) == 2) {

      order2_dist.addValue(d_concentration[j].value,1.0);

    } else {

      order3g_dist.addValue(d_concentration[j].value,1.0);
    }
  }
  for(unsigned int j=0; j<d_temperature.size(); ++j) {
    temp_dist.addValue(d_temperature[j].value,1.0);
  }
  // print out current derivative terms
  CVodeGetCurrentStep(idtCtrl->cvodeCtrl.cvodeMemPtr,&dt_last);
#if defined SUNDIALS2 || defined SUNDIALS3
  CVSpilsGetNumPrecEvals(idtCtrl->cvodeCtrl.cvodeMemPtr,&num_jac);
#else
  CVodeGetNumPrecEvals(idtCtrl->cvodeCtrl.cvodeMemPtr,&num_jac);
#endif
  fprintf(fptr,
          "# gnuplot index %d: (CVSpilsGetNumPrecEvals = %d)\n",
          (int)num_jac-1,
          (int)num_jac);
  fprintf(fptr,
          "# current ODE time    [s]: %.18g\n",
          t_current);
  fprintf(fptr,
          "# last ODE time step  [s]: %.18g\n",
          dt_last);
  fprintf(fptr,
          "# current temperature [K]: %.18g\n",
          thermo_state.temperature);
  fprintf(fptr,
          "# current pressure   [Pa]: %.18g\n",
          thermo_state.pressure);
  fprintf(fptr,
          "# column 1 [Hz]: mininum bin frequency\n");
  fprintf(fptr,
          "# column 2 [Hz]: maximum bin frequency\n");
  fprintf(fptr,
          "# column 3 [#]: d/dC first order steps (no 3rd body or falloff) Total: %16.9e\n",
          order1_dist.total_weight());
  fprintf(fptr,
          "# column 4 [#]: d/dC second order steps (no 3rd body or falloff) Total: %16.9e\n",
          order2_dist.total_weight());
  fprintf(fptr,
          "# column 5 [#]: d/dC third order and higher steps (no 3rd body or falloff) Total: %16.9e\n",
          order3g_dist.total_weight());
  fprintf(fptr,
          "# column 6 [#]: d/dC third-body and falloff steps (all orders) Total: %16.9e\n",
          thirdfall_dist.total_weight());
  fprintf(fptr,
          "# column 7 [#]: d/dT of steps (all orders) Total: %16.9e\n",
          temp_dist.total_weight());
  fprintf(fptr,
          "%16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e\n",
          -1.0e300,
          temp_dist.min_range(),
	  order1_dist.under_range_weight(),
	  order2_dist.under_range_weight(),
	  order3g_dist.under_range_weight(),
	  thirdfall_dist.under_range_weight(),
	  temp_dist.under_range_weight());

  for(int j=0; j<temp_dist.num_bins(); ++j) {
    fprintf(fptr,
            "%16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e\n",
            temp_dist.getBinMin(j),
	    temp_dist.getBinMax(j),
	    order1_dist.getBinWeight(j),
	    order2_dist.getBinWeight(j),
	    order3g_dist.getBinWeight(j),
	    thirdfall_dist.getBinWeight(j),
	    temp_dist.getBinWeight(j));
  }
  fprintf(fptr,
          "%16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e\n",
          temp_dist.max_range(),
          1.0e300,
	  order1_dist.over_range_weight(),
	  order2_dist.over_range_weight(),
	  order3g_dist.over_range_weight(),
	  thirdfall_dist.over_range_weight(),
	  temp_dist.over_range_weight());
  fprintf(fptr,"\n\n"); // gnuplot index record separation


  (*order1_cum)    += order1_dist;
  (*order2_cum)    += order2_dist;
  (*order3g_cum)   += order3g_dist;
  (*thirdfall_cum) += thirdfall_dist;
  (*temp_cum)      += temp_dist;

}
