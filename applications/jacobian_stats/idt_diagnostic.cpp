#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fenv.h>
#ifdef _WIN32
#include <io.h>
#define F_OK 00
#define R_OK 04
#define access _access
#else
#include "unistd.h"
#endif

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


#include "zerork/mechanism.h"
#include "utilities/distribution.h"
#include "transport/binary_collision.h"

#include "event_counter.h"
#include "atol_crossing.h"
#include "BasicReactorIFP.h"
#include "mechanism_stats.h"
#include "jacobian_stats.h"
#include "idtSolvers.h"
#include "idt_diagnostic.h"
#include "rhs_tests.h"

int main(int argc, char *argv[])
{
  //feenableexcept(FE_INVALID | FE_OVERFLOW); // may throw exception
  // or halt program

  int found_idt = 0;
  double delay_temp, delay_time;

  const int num_freq_bins=28;
  checkCommandLine(argc,argv);
  long int num_jac_last, num_jac_current;
  double t_current,t_target,dt_last;
  double cpu_time;
  double *idt_results;

  idtControlParams idtCtrl(argv[1],
                           1); // 1=write the mechanism parser log file
  BasicReactorIFP parser(argv[1]);

  printf("idtControlParams construction complete\n");  fflush(stdout);
  zerork::utilities::Distribution dt_dist(18,
                               1e-20,
                               1e-2,
                               true);
  zerork::utilities::Distribution order1_cum(num_freq_bins,
                               1e-8,
                               1e20,
                               true);
  zerork::utilities::Distribution order2_cum(order1_cum);
  zerork::utilities::Distribution order3g_cum(order1_cum);
  zerork::utilities::Distribution thirdfall_cum(order1_cum);
  zerork::utilities::Distribution temp_cum(order1_cum);

  CVReactorState cv_state;
  FILE *odeProbFile=fopen(idtCtrl.integratorProbFile.c_str(),"w");
  FILE *jacDistFile=fopen(idtCtrl.jacobianStatFile.c_str(),"w");
  FILE *jacRawFile=fopen(idtCtrl.jacobianRawFile.c_str(),"w");
  FILE *outFile=fopen(idtCtrl.outFile.c_str(),"w");

  // negative species stats
  double most_negative = 1.0e300;
  double most_negative_time = 0.0;
  double most_negative_temperature = 0.0;
  int most_negative_id = 0;
  EventCounter negative_species_counter;
  std::vector<int> events, counts;

  // ode problem species stats
  std::vector<double> cvode_species_error;
  std::vector<int> cvode_species_id;
  double max_cvode_species_error;
  int max_cvode_species_id;
  EventCounter cvode_species_counter;
  // ode problem reaction stats
  std::vector<double> cvode_reaction_error;
  std::vector<int> cvode_reaction_id;
  double max_cvode_reaction_error;
  int max_cvode_reaction_id;
  EventCounter cvode_reaction_counter;
  AbsoluteToleranceCrossing atol_crossing(parser.absTol());
  TimeHistoryElement history_element;
  std::vector<TimeHistoryElement> error_history;
  transport::BinaryCollision binary_collision(idtCtrl.mech,
                                              idtCtrl.scan_temperature_ref);

  if(odeProbFile == NULL) {
    printf("ERROR: could not open file %s for write operation.\n",
           idtCtrl.integratorProbFile.c_str());
    exit(-1);
  }
  writeMechanismReport(&parser,&idtCtrl);

  printf("# Initial Composition:\n");
  for(int j=0; j<idtCtrl.mech->getNumSpecies(); ++j) {
    if(idtCtrl.initMoleFrac[j] > 0.0) {
      printf("#   species %4d %16s: %14.7e [by vol], %14.7e [by mass]\n",
             j,
             idtCtrl.mech->getSpeciesName(j),
             idtCtrl.initMoleFrac[j],
             idtCtrl.initMassFrac[j]);
      fflush(stdout);
    }
  }

  idt_results = new double[idtCtrl.cvodeCtrl.nRoots];

  initializeIdtStep(&idtCtrl,
                    idt_results);

  t_current=0.0;
  num_jac_last = num_jac_current = 0;

  int next_step_check = checkNextIdtStep(&idtCtrl, t_current, 1);
  if(next_step_check != 0) {
    fclose(odeProbFile);
    fclose(jacDistFile);
    fclose(jacRawFile);
    fclose(outFile);
    delete [] idt_results;
    printf("ERROR: isfinite() failed for state vector or derivatives on the first step.\n");
    printf("       Exiting now!\n");
    fflush(stdout);
    exit(-1);
  }

  getState(&idtCtrl,&cv_state);
  atol_crossing.UpdateCurrentState(cv_state);
  delay_time = 1.0e300;
  delay_temp = cv_state.temperature;

  // set the step limiter
  if(parser.use_bimolecular_limit()) {

    binary_collision.GetCollisionRateCoefficient(cv_state.temperature,
     		         &idtCtrl.odeUserParams.bimolecular_limiter[0]);
    for(int j=0; j< idtCtrl.nStep; ++j) {
      if(idtCtrl.odeUserParams.bimolecular_limiter[j] > 0.0) {
        idtCtrl.odeUserParams.step_limiter[j] =
          idtCtrl.odeUserParams.bimolecular_limiter[j]*
          parser.bimolecular_limit();
      }
    }
  }
  if(parser.use_unimolecular_limit()) {
    //printf("# INFO: use unimolecular limit = %14.7e\n",
    //      parser.unimolecular_limit());
    // unimolecular should only need to be set once
    for(int j=0; j< idtCtrl.nStep; ++j) {
      if(idtCtrl.odeUserParams.unimolecular_limiter[j] > 0.0) {
        idtCtrl.odeUserParams.step_limiter[j] =
          idtCtrl.odeUserParams.unimolecular_limiter[j];
        //printf("[%4d]: %14.7e unimolecular\n",j,
        //       idtCtrl.odeUserParams.step_limiter[j]);
      }
    }
  }

  //t_next_print = idtCtrl.printDeltaTime;
  while(t_current < idtCtrl.cvodeCtrl.maxTime) {
    t_target = idtCtrl.cvodeCtrl.maxTime;
    int march_flag = marchIdtOneStep(&idtCtrl,
                                  t_current,
                                  &t_target,
                                  &cpu_time);
    CVodeGetCurrentStep(idtCtrl.cvodeCtrl.cvodeMemPtr,&dt_last);
#if defined SUNDIALS2 || defined SUNDIALS3
    CVSpilsGetNumPrecEvals(idtCtrl.cvodeCtrl.cvodeMemPtr,&num_jac_current);
#else
    CVodeGetNumPrecEvals(idtCtrl.cvodeCtrl.cvodeMemPtr,&num_jac_current);
#endif
    if(num_jac_current > num_jac_last) {
      printf("#### NEW JAC = %d ####\n",(int)num_jac_current);
      updateROPDerivativeDistribution(&idtCtrl,
                                      t_current,
                                      &order1_cum,
                                      &order2_cum,
                                      &order3g_cum,
                                      &thirdfall_cum,
                                      &temp_cum,
                                      jacDistFile);
      //getROPDerivativeVector(&idtCtrl,
      //                       &jacobian_terms,
      //                       &min_concentration);
      writeROPDerivativeReport(&idtCtrl,
                               t_current,
                               parser.min_jacobian_raw(),
                               jacRawFile);
      writeROPTempDerivativeReport(&idtCtrl,
                               t_current,
                               parser.min_df_dtemp_raw(),
                               jacRawFile);
    }

    //dt_last = t_target - t_current;
    dt_dist.addValue(dt_last,1.0);
    t_current = t_target;
    num_jac_last = num_jac_current;
    checkAllFlagsFromIdtStep(&idtCtrl,
                             t_current,
                             march_flag,
                             idt_results,
                             odeProbFile);
    getState(&idtCtrl,&cv_state);
    atol_crossing.UpdateCurrentState(cv_state);
    // set the step limiter
    if(parser.use_bimolecular_limit()) {

      binary_collision.GetCollisionRateCoefficient(cv_state.temperature,
       		         &idtCtrl.odeUserParams.bimolecular_limiter[0]);
      for(int j=0; j< idtCtrl.nStep; ++j) {
        if(idtCtrl.odeUserParams.bimolecular_limiter[j] > 0.0) {
          idtCtrl.odeUserParams.step_limiter[j] =
            idtCtrl.odeUserParams.bimolecular_limiter[j]*
            parser.bimolecular_limit();
        }
      }
    }

    history_element.current_time_ = t_current;
    history_element.current_temperature_ = cv_state.temperature;
    history_element.last_time_step_ = dt_last;

    int min_species_id;
    double min_mole_fraction = GetMinMoleFraction(cv_state,
                                                  &min_species_id);
    double negative_mole_fraction_sum = GetNegativeMoleFractionSum(cv_state);

    if(min_mole_fraction < most_negative) {
      most_negative    = min_mole_fraction;
      most_negative_id = min_species_id;
      most_negative_time = t_current;
      most_negative_temperature = cv_state.temperature;
    }
    history_element.most_negative_species_id_ = most_negative_id;

    if(min_mole_fraction < 0.0) {
      negative_species_counter.AddEvent(min_species_id);
    }
    // Check cvode scaled errors by species
    history_element.error_species_id_.clear();
    GetCVodeErrorSpecies(&idtCtrl,
                         parser.min_cvode_species_error(),
                         &cvode_species_id,
                         &cvode_species_error);
    max_cvode_species_error = -1.0e300;
    for(size_t j=0; j<cvode_species_id.size(); ++j) {
      if(max_cvode_species_error < cvode_species_error[j]) {
        max_cvode_species_error = cvode_species_error[j];
        max_cvode_species_id = cvode_species_id[j];
      }
      if(cvode_species_error[j] >= parser.min_cvode_species_error()) {
        cvode_species_counter.AddEvent(cvode_species_id[j]);
        history_element.error_species_id_.push_back(cvode_species_id[j]);
      }
    }

    // Check cvode scaled errors by reactions
    GetCVodeErrorReactions(&idtCtrl,
                         parser.min_cvode_reaction_error(),
                         &cvode_reaction_id,
                         &cvode_reaction_error);
    max_cvode_reaction_error = -1.0e300;
    for(size_t j=0; j<cvode_reaction_id.size(); ++j) {
      if(max_cvode_reaction_error < cvode_reaction_error[j]) {
        max_cvode_reaction_error = cvode_reaction_error[j];
        max_cvode_reaction_id = cvode_reaction_id[j];
      }
      if(cvode_reaction_error[j] >= parser.min_cvode_reaction_error()) {
        cvode_reaction_counter.AddEvent(cvode_reaction_id[j]);
      }
    }

    error_history.push_back(history_element);

    printf("%14.7e  %10.3e  %10.4f  %10.3e  %10.3e  %5d  %-16s  %10.3e  %5d  %-16s  %10.3e  %s\n",
           t_current,
           dt_last,
           cv_state.temperature,
           negative_mole_fraction_sum,
           min_mole_fraction,
           min_species_id,
           idtCtrl.mech->getSpeciesName(min_species_id),
	   max_cvode_species_error,
           max_cvode_species_id,
           idtCtrl.mech->getSpeciesName(max_cvode_species_id),
           max_cvode_reaction_error,
           idtCtrl.mech->getReactionName(max_cvode_reaction_id));
    fprintf(outFile,
            "%14.7e  %10.3e  %10.4f  %10.3e  %10.3e  %5d  %-16s  %10.3e  %5d  %-16s  %10.3e  %s\n",
            t_current,
            dt_last,
            cv_state.temperature,
            negative_mole_fraction_sum,
            min_mole_fraction,
            min_species_id,
	    idtCtrl.mech->getSpeciesName(min_species_id),
	    max_cvode_species_error,
            max_cvode_species_id,
            idtCtrl.mech->getSpeciesName(max_cvode_species_id),
            max_cvode_reaction_error,
            idtCtrl.mech->getReactionName(max_cvode_reaction_id));

     if(cv_state.temperature >= idtCtrl.initTemp+400.0 && found_idt == 0) {
      found_idt = 1;
      delay_time = t_current;
      delay_temp = cv_state.temperature;
    }
    fflush(stdout);
  }

  printf("\n\n# Time step distribution:\n");
  printf("# Total Number of Steps      = %.18g\n",dt_dist.total_weight());
  printf("# Total lower than minimum   = %.18g\n",dt_dist.under_range_weight());
  printf("# Total greater than maximum = %.18g\n",dt_dist.over_range_weight());
  for(int j=0; j<dt_dist.num_bins(); ++j) {
    printf("%.18g  %.18g  %.18g  %.18g\n",
           dt_dist.getBinMin(j),
           dt_dist.getBinMax(j),
           dt_dist.getBinWeight(j),
           dt_dist.getBinFraction(j));
      }

  fprintf(outFile,"\n\n# Time step distribution:\n");
  fprintf(outFile,"# Total Number of Steps      = %.18g\n",dt_dist.total_weight());
  fprintf(outFile,"# Total lower than minimum   = %.18g\n",dt_dist.under_range_weight());
  fprintf(outFile,"# Total greater than maximum = %.18g\n",dt_dist.over_range_weight());
  for(int j=0; j<dt_dist.num_bins(); ++j) {
    fprintf(outFile,"%.18g  %.18g  %.18g  %.18g\n",
           dt_dist.getBinMin(j),
           dt_dist.getBinMax(j),
           dt_dist.getBinWeight(j),
           dt_dist.getBinFraction(j));
      }
  // print out the cumulative Jacobian term distribution
  fprintf(jacDistFile,
          "# gnuplot index %d: (cumulative over entire idt calculation\n",
          (int)num_jac_current);
  fprintf(jacDistFile,
          "# column 1 [Hz]: mininum bin frequency\n");
  fprintf(jacDistFile,
          "# column 2 [Hz]: maximum bin frequency\n");
   fprintf(jacDistFile,
          "# column 3 [#]: d/dC first order steps (no 3rd body or falloff) Total: %16.9e\n",
          order1_cum.total_weight());
  fprintf(jacDistFile,
          "# column 4 [#]: d/dC second order steps (no 3rd body or falloff) Total: %16.9e\n",
          order2_cum.total_weight());
  fprintf(jacDistFile,
          "# column 5 [#]: d/dC third order and higher steps (no 3rd body or falloff) Total: %16.9e\n",
          order3g_cum.total_weight());
  fprintf(jacDistFile,
          "# column 6 [#]: d/dC third-body and falloff steps (all orders) Total: %16.9e\n",
          thirdfall_cum.total_weight());
  fprintf(jacDistFile,
          "# column 7 [#]: d/dT of steps (all orders) Total: %16.9e\n",
          temp_cum.total_weight());
  fprintf(jacDistFile,
          "%16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e\n",
          -1.0e300,
          temp_cum.min_range(),
	  order1_cum.under_range_weight(),
	  order2_cum.under_range_weight(),
	  order3g_cum.under_range_weight(),
	  thirdfall_cum.under_range_weight(),
	  temp_cum.under_range_weight());

  for(int j=0; j<temp_cum.num_bins(); ++j) {
    fprintf(jacDistFile,
            "%16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e\n",
            temp_cum.getBinMin(j),
	    temp_cum.getBinMax(j),
	    order1_cum.getBinWeight(j),
	    order2_cum.getBinWeight(j),
	    order3g_cum.getBinWeight(j),
	    thirdfall_cum.getBinWeight(j),
	    temp_cum.getBinWeight(j));
  }
  fprintf(jacDistFile,
          "%16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e\n",
          temp_cum.max_range(),
          1.0e300,
	  order1_cum.over_range_weight(),
	  order2_cum.over_range_weight(),
	  order3g_cum.over_range_weight(),
	  thirdfall_cum.over_range_weight(),
	  temp_cum.over_range_weight());

  printf("\n\n# Most negative mole fraction = %.18g\n",
         most_negative);
  printf("#    species id     = %d\n",
         most_negative_id);
  printf("#    species name   = %s\n",
         idtCtrl.mech->getSpeciesName(most_negative_id));
  printf("#    at time        = %.18g\n",
         most_negative_time);
  printf("#    at temperature = %.18g\n",
         most_negative_temperature);
  fprintf(outFile,"\n\n# Most negative mole fraction = %.18g\n",
          most_negative);
  fprintf(outFile,"#    species id     = %d\n",
          most_negative_id);
  fprintf(outFile,"#    species name   = %s\n",
          idtCtrl.mech->getSpeciesName(most_negative_id));
  fprintf(outFile,"#    at time        = %.18g\n",
          most_negative_time);
  fprintf(outFile,"#    at temperature = %.18g\n",
          most_negative_temperature);

  negative_species_counter.GetSortedEventList(&events,
                                              &counts);
  printf("# Cumulative statistics for the most negative mole fractions\n"
         "# Column 1: species index\n"
         "# Column 2: species name\n"
         "# Column 3: number of times the species was the most negative mole fraction\n");
  fprintf(outFile,
    "# Cumulative statistics for the most negative mole fractions"
    "# Column 1: species index\n"
    "# Column 2: species name\n"
    "# Column 3: number of times the species was the most negative mole fraction\n");
  for(int j=0; j<negative_species_counter.GetNumDistinctEvents(); ++j) {
    printf("%5d  %16s  %6d\n",
           events[j],
           idtCtrl.mech->getSpeciesName(events[j]),
           counts[j]);
    fprintf(outFile,
            "%5d  %16s  %6d\n",
            events[j],
            idtCtrl.mech->getSpeciesName(events[j]),
            counts[j]);
  }

  cvode_species_counter.GetSortedEventList(&events,
                                           &counts);
  printf("\n\n# Cumulative statistics for the cvode error species. This is\n"
         "# the number of the times the cvode scaled error measure\n"
         "# was equal to or exceeded %.18g. Note that a measure of one\n"
         "# is deemed to just satisfy the relative and absolute error\n"
         "# tolerances.\n"
         "#     Column  1: species index\n"
         "#     Column  2: species name\n"
         "#     Column  3: number of times the species exceed the scaled error threshold\n"
         "#     Column  4: [s] min time of cvode error\n"
         "#     Column  5: [s] max time of cvode error\n"
         "#     Column  6: [K] min temperature of cvode error\n"
         "#     Column  7: [K] max temperature of cvode error\n"
         "#     Column  8: [s] min time step of cvode error\n"
         "#     Column  9: [s] max time step of cvode error\n"
         "#     Column 10: [s] mean time step of cvode error\n",
         parser.min_cvode_species_error());
  fprintf(outFile,
          "\n\n# Cumulative statistics for the cvode error species. This is\n"
         "# the number of the times the cvode scaled error measure\n"
         "# was equal to or exceeded %.18g. Note that a measure of one\n"
         "# is deemed to just satisfy the relative and absolute error\n"
         "# tolerances.\n"
         "#     Column  1: species index\n"
         "#     Column  2: species name\n"
         "#     Column  3: number of times the species exceed the scaled error threshold\n"
         "#     Column  4: [s] min time of cvode error\n"
         "#     Column  5: [s] max time of cvode error\n"
         "#     Column  6: [K] min temperature of cvode error\n"
         "#     Column  7: [K] max temperature of cvode error\n"
         "#     Column  8: [s] min time step of cvode error\n"
         "#     Column  9: [s] max time step of cvode error\n"
	  "#     Column 10: [s] exp(mean log(time step)) of cvode error\n",
	  parser.min_cvode_species_error());
  for(int j=0; j<cvode_species_counter.GetNumDistinctEvents(); ++j) {
    double e_min_time, e_max_time, e_min_temperature, e_max_temperature;
    double e_min_time_step, e_max_time_step, e_mean_time_step;
    int num_errors = GetErrorHistoryRanges(error_history,
                                           events[j],
                                           &e_min_time,
                                           &e_max_time,
                                           &e_min_temperature,
                                           &e_max_temperature,
                                           &e_min_time_step,
                                           &e_max_time_step,
                                           &e_mean_time_step);
    if(num_errors != counts[j]) {
      printf("# ERROR: event count does not match the error history list\n");
      printf("         for species id = %d (%s)\n",
             events[j],
             idtCtrl.mech->getSpeciesName(events[j]));
    }

    printf("%5d  %-16s  %6d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e\n",
           events[j],
           idtCtrl.mech->getSpeciesName(events[j]),
           counts[j],
           e_min_time,
           e_max_time,
           e_min_temperature,
           e_max_temperature,
           e_min_time_step,
           e_max_time_step,
           e_mean_time_step);
    fprintf(outFile,
            "%5d  %-16s  %6d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e\n",
            events[j],
            idtCtrl.mech->getSpeciesName(events[j]),
            counts[j],
            e_min_time,
            e_max_time,
            e_min_temperature,
            e_max_temperature,
            e_min_time_step,
            e_max_time_step,
            e_mean_time_step);
  }

  atol_crossing.GetSortedCrossingList(&events,
                                      &counts);
  printf("\n\n# Cumulative statistics for the absolute value of the species\n"
         "# mass fractions crossing the absolute tolerance = %10.3e three or more times\n"
         "# Column 1: species index\n"
         "# Column 2: species name\n"
         "# Column 3: number of times the absolute value of the species mass\n"
         "#           fraction crossed the absolute tolerance\n",
         parser.absTol());
  fprintf(outFile,
          "\n\n# Cumulative statistics for the absolute value of the species\n"
          "# mass fractions crossing the absolute tolerance = %10.3e three or more times\n"
          "# Column 1: species index\n"
          "# Column 2: species name\n"
          "# Column 3: number of times the absolute value of the species mass\n"
          "#           fraction crossed the absolute tolerance\n",
          parser.absTol());

  for(int j=0; j<(int)events.size(); ++j) {
    if(counts[j] > 2) {
      printf("%5d  %16s  %6d\n",
             events[j],
             idtCtrl.mech->getSpeciesName(events[j]),
             counts[j]);
      fprintf(outFile,
              "%5d  %16s  %6d\n",
              events[j],
              idtCtrl.mech->getSpeciesName(events[j]),
              counts[j]);
    }
  }


  // cvode_reaction_counter.GetSortedEventList(&events,
  //                                           &counts);
  // printf("\n\n# Cumulative statistics for the cvode error by reaction. This is\n"
  //        "# the number of the times the cvode scaled error measure\n"
  //        "# was equal to or exceeded %.18g. Note that a measure of one\n"
  //        "# is deemed to just satisfy the relative and absolute error\n"
  //        "# tolerances.\n"
  //        "#     Column 1: reaction index\n"
  //        "#     Column 2: reaction name\n"
  //        "#     Column 3: number of times the species exceed the scaled error threshold\n",
  //        parser.min_cvode_reaction_error());
  // fprintf(outFile,
  //         "\n\n# Cumulative statistics for the cvode error by reaction. This is\n"
  //        "# the number of the times the cvode scaled error measure\n"
  //        "# was equal to or exceeded %.18g. Note that a measure of one\n"
  //        "# is deemed to just satisfy the relative and absolute error\n"
  //        "# tolerances.\n"
  //        "#     Column 1: reaction index\n"
  //        "#     Column 2: reaction name\n"
  //        "#     Column 3: number of times the species exceed the scaled error threshold\n",
  // 	  parser.min_cvode_reaction_error());
  // for(int j=0; j<cvode_reaction_counter.GetNumDistinctEvents(); ++j) {
  //   printf("%5d  %-16s  %6d\n",
  //          events[j],
  //          idtCtrl.mech->getReactionName(events[j]),
  //          counts[j]);
  //   fprintf(outFile,
  //           "%5d  %-16s  %6d\n",
  //           events[j],
  //           idtCtrl.mech->getReactionName(events[j]),
  //           counts[j]);
  // }


  printf("\n\n# Initial Temperature [K]: %9.4f\n",idtCtrl.initTemp);
  printf("# Delay   Temperature [K]: %9.4f\n",delay_temp);
  printf("# Delay   Time        [s]: %16.7e\n",delay_time);
  fprintf(outFile,"\n\n# Initial Temperature [K]: %9.4f\n",idtCtrl.initTemp);
  fprintf(outFile,"# Delay   Temperature [K]: %9.4f\n",delay_temp);
  fprintf(outFile,"# Delay   Time        [s]: %16.7e\n",delay_time);
  fclose(odeProbFile);
  fclose(jacDistFile);
  fclose(jacRawFile);
  fclose(outFile);
  delete [] idt_results;
  return 0;
}

void getInitialActiveSpecies(BasicReactorIFP *parser,
                             zerork::mechanism *mech,
                             std::vector<int> *initial_active_species)
{
  int species_id;
  typedef map< string, double > comp_t;
  comp_t::const_iterator iter;

  std::vector<int> active_species;
  active_species.assign(mech->getNumSpecies(),0);

  for(iter =  parser->fuelComp().begin();
      iter != parser->fuelComp().end();
      ++iter) {

    species_id = mech->getIdxFromName(iter->first.c_str());
    if(species_id == -1) { 
       printf("ERROR: fuel species %s: not found in mechanism.\n", iter->first.c_str());
       exit(1);
    }
    if(iter->second > 0.0) {
      active_species[species_id]=1;
    }
  }
  for(iter =  parser->oxidizerComp().begin();
      iter != parser->oxidizerComp().end();
      ++iter) {

    species_id = mech->getIdxFromName(iter->first.c_str());
    if(species_id == -1) { 
       printf("ERROR: oxid species %s: not found in mechanism.\n", iter->first.c_str());
       exit(1);
    }
    if(iter->second > 0.0) {
      active_species[species_id]=1;
    }
  }

  initial_active_species->clear();
  for(int j=0; j<mech->getNumSpecies(); ++j) {
    if(active_species[j] == 1) {
      initial_active_species->push_back(j);
    }
  }
}

void checkCommandLine(int inpArgc, char *inpArgv[])
{
  if(inpArgc != 2) {
    printf("ERROR: incorrect command line usage\n");
    printf("       use instead %s <input file>\n",
           inpArgv[0]);
    exit(-1);
  }
  if(access(inpArgv[1],R_OK) != 0) {
    printf("ERROR: can not open input file %s for read access\n",
           inpArgv[1]);
    exit(-2);
  }
}

void writeMechanismReport(BasicReactorIFP *parser,
                          idtControlParams *idt_control)
{
  std::string info;
  std::vector<int> initial_active_species;
  FILE *fptr=fopen(idt_control->mechStatFile.c_str(),"w");

  if(fptr == NULL) {
    printf("ERROR: In writeMechanismReport(...),\n");
    printf("       could not open file %s for write operation.\n",
           idt_control->mechStatFile.c_str());
    fflush(stdout);
    exit(-1);
  }

  fprintfDashedLine(fptr);
  fprintf(fptr,"# Mechanism file     : %s\n",idt_control->mechFile.c_str());
  fprintf(fptr,"# Thermodynamics file: %s\n",idt_control->thermFile.c_str());
  fprintf(fptr,"# Conversion log file: %s\n",idt_control->mechLogFile.c_str());

  fprintfDashedLine(fptr);
  getBasicMechanismReport(idt_control->mech,&info);
  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  getSpeciesSummaryReport(idt_control->mech,&info);
  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  getInertSpeciesReport(idt_control->mech,&info);
  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  getAdductSpeciesReport(idt_control->mech,&info);
  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  getLoneReactionSpeciesReport(idt_control->mech,&info);
  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  for(int j=0; j<idt_control->mech->getNumSpecies(); ++j) {
    if(idt_control->initMoleFrac[j] > 0.0) {
      initial_active_species.push_back(j);
    }
  }
//  getInactiveReport(&initial_active_species,
//                    idt_control->mech,
//                    &info);
//  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  getSourceSpeciesReport(idt_control->mech,&info);
  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  getSinkSpeciesReport(idt_control->mech,&info);
  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  getNegativeRateReport(idt_control->mech,
                        idt_control->initPres,
                        idt_control->min_scan_temperature,
                        idt_control->max_scan_temperature,
                        idt_control->num_scans,
                        &info);
  fprintf(fptr,"%s",info.c_str());

  fprintfDashedLine(fptr);
  getMaxUnimolecularRateReport(idt_control->mech,
                               idt_control->initPres,
                               idt_control->min_scan_temperature,
                               idt_control->max_scan_temperature,
                               idt_control->num_scans,
                               &info);
  fprintf(fptr,"%s",info.c_str());
  fprintfDashedLine(fptr);
  getSortedUnimolecularRateReport(idt_control->mech,
                                  idt_control->initPres, // reference pressure
                                  idt_control->min_scan_temperature,
                                  idt_control->max_scan_temperature,
                                  idt_control->num_scans,
                                  parser->min_unimolecular_rate(),
                                  &info);
  fprintf(fptr,"%s",info.c_str());

  transport::BinaryCollision binary_collision_chk(idt_control->mech,
                                            idt_control->scan_temperature_ref);
  fprintfDashedLine(fptr);
  binary_collision_chk.GetSpeciesReport(&info);
  fprintf(fptr,"%s",info.c_str());
  fprintfDashedLine(fptr);
  binary_collision_chk.GetStepProbabilityReport(idt_control->initPres,
                                             idt_control->min_scan_temperature,
                                             idt_control->max_scan_temperature,
					     idt_control->num_scans,
					     parser->min_bimolecular_prob(),
                                             &info);
  fprintf(fptr,"%s",info.c_str());


  fclose(fptr);
}

// return the number of species found above the error threshold
int GetCVodeErrorSpecies(idtControlParams *idt_control,
                         const double threshold,
                         std::vector<int> *species_id,
                         std::vector<double> *species_error)
{
  const int num_species = idt_control->nSpc;
  const int num_states  = num_species+2;

  N_Vector error_weights;
  N_Vector local_errors;

  error_weights = N_VNew_Serial(num_states);
  local_errors  = N_VNew_Serial(num_states);

  double max_error = -1.0e300;
  int max_species_id = -1;
  int num_over_threshold;

  CVodeGetErrWeights(idt_control->cvodeCtrl.cvodeMemPtr,
                     error_weights);

  CVodeGetEstLocalErrors(idt_control->cvodeCtrl.cvodeMemPtr,
                         local_errors);
  species_id->clear();
  species_error->clear();

  for(int j=0; j<num_species; ++j) {
    double scaled_error = fabs(NV_Ith_S(error_weights,j)*
                               NV_Ith_S(local_errors,j));
    if(scaled_error > max_error) {
      max_error = scaled_error;
      max_species_id = j;
    }
    if(scaled_error >= threshold) {
      species_id->push_back(j);
      species_error->push_back(scaled_error);
    }
  }
  num_over_threshold = (int)species_id->size();
  if(num_over_threshold == 0) {
    // store the species with the max error
    species_id->push_back(max_species_id);
    species_error->push_back(max_error);
  }
  N_VDestroy_Serial(error_weights);
  N_VDestroy_Serial(local_errors);
  return num_over_threshold;
}
// return the number of reactions found above the error threshold
int GetCVodeErrorReactions(idtControlParams *idt_control,
                           const double threshold,
                           std::vector<int> *reaction_id,
                           std::vector<double> *reaction_error)
{
  const int num_species   = idt_control->nSpc;
  const int num_states    = num_species+2;
  const int num_reactions = idt_control->nRxn;

  std::vector<double> species_error(num_species);
  N_Vector error_weights;
  N_Vector local_errors;

  error_weights = N_VNew_Serial(num_states);
  local_errors  = N_VNew_Serial(num_states);

  double max_error = -1.0e300;
  int max_reaction_id = -1;
  int num_over_threshold;

  CVodeGetErrWeights(idt_control->cvodeCtrl.cvodeMemPtr,
                     error_weights);

  CVodeGetEstLocalErrors(idt_control->cvodeCtrl.cvodeMemPtr,
                         local_errors);
  reaction_id->clear();
  reaction_error->clear();

  for(int j=0; j<num_species; ++j) {
    species_error[j] = fabs(NV_Ith_S(error_weights,j)*
                            NV_Ith_S(local_errors,j));
  }
  for(int j=0; j<num_reactions; ++j) {
    double reaction_error_sum = 0.0;
    int step_index = idt_control->mech->getStepIdxOfRxn(j, 1); // forward step
    int num_reactants = idt_control->mech->getOrderOfStep(step_index);
    int num_products  = idt_control->mech->getNumProductsOfStep(step_index);
    int participant_id;
    int num_error_species = 0;
    for(int k=0; k<num_reactants; ++k) {
      participant_id = idt_control->mech->getSpecIdxOfStepReactant(step_index,
                                                                   k);
      reaction_error_sum += species_error[participant_id];
      if(threshold < species_error[participant_id]) {
        ++num_error_species;
      }
    }
    for(int k=0; k<num_products; ++k) {
      participant_id = idt_control->mech->getSpecIdxOfStepProduct(step_index,
                                                                  k);
      reaction_error_sum += species_error[participant_id];
      if(threshold < species_error[participant_id]) {
        ++num_error_species;
      }
    }
    if(num_error_species < 2) {
      reaction_error_sum = 0.0;
    }

    if(max_error < reaction_error_sum) {
      max_error = reaction_error_sum;
      max_reaction_id = j;
    }
    if(threshold < reaction_error_sum) {
      reaction_id->push_back(j);
      reaction_error->push_back(reaction_error_sum);
    }
  }

  num_over_threshold = (int)reaction_id->size();
  if(num_over_threshold == 0) {
    // store the species with the max error
    reaction_id->push_back(max_reaction_id);
    reaction_error->push_back(max_error);
  }
  N_VDestroy_Serial(error_weights);
  N_VDestroy_Serial(local_errors);
  return num_over_threshold;
}

int GetErrorHistoryRanges(const std::vector<TimeHistoryElement> &error_history,
                          const int species_id,
                          double *min_time,
                          double *max_time,
                          double *min_temperature,
                          double *max_temperature,
                          double *min_time_step,
                          double *max_time_step,
                          double *mean_time_step)
{
  int num_species_errors = 0;
  *min_time = 1.0e300;
  *max_time =-1.0e300;
  *min_temperature = 1.0e300;
  *max_temperature =-1.0e300;
  *min_time_step = 1.0e300;
  *max_time_step =-1.0e300;
  *mean_time_step = 0.0;

  for(size_t j=0; j<error_history.size(); ++j) {
    for(size_t k=0; k<error_history[j].error_species_id_.size(); ++k) {
      if(species_id == error_history[j].error_species_id_[k]) {

        if(error_history[j].current_time_ < *min_time) {
          *min_time = error_history[j].current_time_;
        }
        if(error_history[j].current_time_ > *max_time) {
          *max_time = error_history[j].current_time_;
        }
        if(error_history[j].current_temperature_ < *min_temperature) {
          *min_temperature = error_history[j].current_temperature_;
        }
        if(error_history[j].current_temperature_ > *max_temperature) {
          *max_temperature = error_history[j].current_temperature_;
        }
        if(error_history[j].last_time_step_ < *min_time_step) {
          *min_time_step = error_history[j].last_time_step_;
        }
        if(error_history[j].last_time_step_ > *max_time_step) {
          *max_time_step = error_history[j].last_time_step_;
        }

        *mean_time_step += log(error_history[j].last_time_step_);
	++num_species_errors;
      }
    }
  }
  if(num_species_errors > 1) {
    *mean_time_step /= (double)num_species_errors;
    *mean_time_step = exp(*mean_time_step);
  }

  return num_species_errors;
}
