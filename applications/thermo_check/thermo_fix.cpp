#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//#include <lapacke.h>

#include "lapack_wrapper.h"
#include "janaf_thermo.h"
#include "thermo_fix.h"

double SpecificHeatDelta(const double T,
                         const double coef_low[],
                         const double coef_high[])
{
  return SpecificHeat(T,coef_high)-SpecificHeat(T,coef_low);
}

double SpecificHeatDerivativeDelta(const double T,
                                   const double coef_low[],
                                   const double coef_high[])
{
  return SpecificHeatDerivative(T,coef_high)-
         SpecificHeatDerivative(T,coef_low);
}

double EnthalpyDelta(const double T,
                     const double coef_low[],
                     const double coef_high[])
{
  return Enthalpy(T,coef_high)-Enthalpy(T,coef_low);
}


double EntropyDelta(const double T,
                    const double coef_low[],
                    const double coef_high[])
{
  return Entropy(T,coef_high)-Entropy(T,coef_low);
}

void RefitKeepHigh(const int num_points,
                   const int T_resolution,
                   const double T_fixed,
                   const double T_min,
                   const double T_match,
                   const double T_max,
                   const double coef_low[],
                   const double coef_high[],
                   double *refit_T_match,
                   double refit_coef_low[],
                   double refit_coef_high[])
{
  bool refit_success = true;
  double T_tol;
  double T_point_constraint[2];
  double Cp_point_constraint[2];
  double *T_orig, *Cp_orig;
  double T_slope_constraint,dCp_slope_constraint;
  double H_fixed, S_fixed;
  double new_T_match=RoundToTenPower(-T_resolution,T_match);
  
  if(refit_T_match != NULL) {
    //printf("refit_T_match != NULL => [1] *refit_T_match = %.18g\n",
    //       *refit_T_match);
    T_tol = pow(10.0,-T_resolution-2.0);
    *refit_T_match = MinSpecificHeatDelta(T_tol,
                                          T_min,
                                          T_match,
                                          T_max,
                                          &coef_low[0],
                                          &coef_high[0]);

    if(*refit_T_match <= T_min || *refit_T_match >= T_max) {
      // don't change the T_match value
      printf("WARNING: attempt to find the best T_match failed\n");
      printf("         keeping original T_match = %.18g\n",T_match); 
      *refit_T_match = T_match;
    }
    //printf("[2] refit_T_match = %.18g\n",*refit_T_match);
    *refit_T_match = RoundToTenPower(-T_resolution,*refit_T_match);
    //printf("[3] refit_T_match = %.18g\n",*refit_T_match);

    new_T_match = *refit_T_match;
  }
  //printf("[4] refit_T_match = %.18g\n",*refit_T_match);
  //printf("new_T_match   = %.18g\n",new_T_match);
  
  // set the constraints on Cp
  T_point_constraint[0] = T_fixed;
  T_point_constraint[1] = new_T_match;
  T_slope_constraint    = new_T_match;

  // use the original form to compute the fixed points
  // TODO: decide if the Cp constraint should be evaluated always on the 
  //       low branch
  Cp_point_constraint[0] = SpecificHeat(T_fixed,
                                        &coef_low[0]);

  // evaluate the high temperature branch at the new matching point 
  Cp_point_constraint[1] = SpecificHeat(new_T_match,
                                        &coef_high[0]);
  dCp_slope_constraint   = SpecificHeatDerivative(new_T_match,
                                                  &coef_high[0]);

  // generate the data to which the new lower branch is fit
  T_orig  = new double[num_points];
  Cp_orig = new double[num_points];

  for(int j=0; j<num_points; ++j) {
    T_orig[j]  = T_min + (new_T_match-T_min)*
      static_cast<double>(j)/(num_points-1.0);
    Cp_orig[j] = SpecificHeat(T_orig[j],&coef_low[0]); 
  }

  // evaluate the new low temperature fits
  refit_success = PolynomialFit(num_points, // number of points generated to fit
                               4,          // polynomial degree
                               2,          // number of fixed point constraints
                               1,          // number of slope constraints
                               T_orig,
                               Cp_orig,
                               T_point_constraint,
                               Cp_point_constraint,
                               &T_slope_constraint,
                               &dCp_slope_constraint,
                               &refit_coef_low[0]);
  if(!refit_success) {
    printf("WARNING: PolynomialFit(...) failed,\n");
    printf("         retaining low temperature coefficients for Cp/R\n");
    for(int j=0; j<5; ++j) {
      refit_coef_low[j] = coef_low[j];
    }
  }
  // copy the original high coefficients
  for(int j=0; j<7; ++j) {
    refit_coef_high[j] = coef_high[j];
  }
  // copy the low coefficients for enthalpy and entropy
  refit_coef_low[5] = coef_low[5];
  refit_coef_low[6] = coef_low[6];

  // calculate the fixed temperature enthalpy and entropy from the original
  // thermo definition
  if(T_fixed < T_match) {
    H_fixed = Enthalpy(T_fixed,&coef_low[0]);
    S_fixed = Entropy( T_fixed,&coef_low[0]);
  } else {
    H_fixed = Enthalpy(T_fixed,&coef_high[0]);
    S_fixed = Entropy( T_fixed,&coef_high[0]);
  }
  // adjust the new thermo definition
  if(T_fixed < new_T_match) {
    // 1. adjust the low temperature coefficient branch
    MatchEnthalpy(T_fixed,H_fixed,&refit_coef_low[0]);
    MatchEntropy( T_fixed,S_fixed,&refit_coef_low[0]);
    // 2. adjust the high temprature coefficient branch to match
    H_fixed = Enthalpy(new_T_match,&refit_coef_low[0]);
    S_fixed = Entropy( new_T_match,&refit_coef_low[0]);

    MatchEnthalpy(new_T_match,H_fixed,&refit_coef_high[0]);
    MatchEntropy( new_T_match,S_fixed,&refit_coef_high[0]);
  } else {
    // 1. adjust the high temperature coefficient branch
    MatchEnthalpy(T_fixed,H_fixed,&refit_coef_high[0]);
    MatchEntropy( T_fixed,S_fixed,&refit_coef_high[0]);
    // 2. adjust the low temprature coefficient branch to match
    H_fixed = Enthalpy(new_T_match,&refit_coef_high[0]);
    S_fixed = Entropy( new_T_match,&refit_coef_high[0]);

    MatchEnthalpy(new_T_match,H_fixed,&refit_coef_low[0]);
    MatchEntropy( new_T_match,S_fixed,&refit_coef_low[0]);
  }

  delete [] T_orig;
  delete [] Cp_orig;
}

void RefitKeepHigh(const int num_points,
                   const int T_resolution,
                   const double T_fixed,
                   const JanafThermoData *original,
                   const bool find_new_match,
		   JanafThermoData *refit)
{
  // copy the temperatures to the new JANAF fit
  refit->T_min   = original->T_min;
  refit->T_match = RoundToTenPower(-T_resolution,original->T_match);
  refit->T_max   = original->T_max;
  // initialize coefficients to zero
  for(int j=0; j<7; ++j) {
    refit->low_coef[j] = 0.0;
    refit->high_coef[j] = 0.0;
  }
  // check if only one temperature range is given
  if(refit->T_match <= refit->T_min) {
    // only in the high range, copy high range coefficients to both ranges
    for(int j=0; j<7; ++j) {
      refit->low_coef[j]  = original->high_coef[j];
      refit->high_coef[j] = original->high_coef[j];
    }
    return;
  } else if(refit->T_match >= refit->T_max) {
    // only in the low range, copy low range coefficients to both ranges
    for(int j=0; j<7; ++j) {
      refit->low_coef[j]  = original->low_coef[j];
      refit->high_coef[j] = original->low_coef[j];
    }
    return;
  }    
 

  //printf("num_points = %d\n",num_points);
  //printf("T_resolution = %d\n",T_resolution);
  //printf("T_fixed      = %.18g\n",T_fixed);
  //printf("find_new_match = %d\n",find_new_match);
  //PrintJanaf(*original);
  //printf("orig %.18g refit %.18g\n",original->T_match,refit->T_match);
  //fflush(stdout);
  if(find_new_match) {

    RefitKeepHigh(num_points,
                  T_resolution, // decimal places for the match temperature
                  T_fixed,
                  original->T_min,
                  original->T_match,
                  original->T_max,
                  original->low_coef,
                  original->high_coef,
                  &refit->T_match,
                  refit->low_coef,
                  refit->high_coef);
  } else {

    RefitKeepHigh(num_points,
                  T_resolution, // decimal places for the match temperature
                  T_fixed,
                  original->T_min,
                  original->T_match,
                  original->T_max,
                  original->low_coef,
                  original->high_coef,
                  NULL,
                  refit->low_coef,
                  refit->high_coef);
    // initial set upon function entry
    //refit->T_match = RoundToTenPower(-T_resolution,original->T_match); 
  }
  

}
void RefitKeepHighGlobalTMatch(const int num_points,
                               const int T_resolution,
                               const double T_fixed,
                               const double T_min,
                               const double T_match,
                               const double T_max,
                               const double coef_low[],
                               const double coef_high[],
                               const double global_T_match,
                               double refit_coef_low[],
                               double refit_coef_high[])
{
  bool refit_success = true;
  double T_point_constraint[2];
  double Cp_point_constraint[2];
  double *T_orig, *Cp_orig;
  double T_slope_constraint,dCp_slope_constraint;
  double H_fixed, S_fixed;
  double new_T_match=RoundToTenPower(-T_resolution,global_T_match);
  //printf("new_T_match   = %.18g\n",new_T_match);
  
  // set the constraints on Cp
  T_point_constraint[0] = T_fixed;
  T_point_constraint[1] = new_T_match;
  T_slope_constraint    = new_T_match;

  // use the original form to compute the fixed points
  // TODO: decide if the Cp constraint should be evaluated always on the 
  //       low branch
  Cp_point_constraint[0] = SpecificHeat(T_fixed,
                                        &coef_low[0]);

  // evaluate the high temperature branch at the new matching point 
  Cp_point_constraint[1] = SpecificHeat(new_T_match,
                                        &coef_high[0]);
  dCp_slope_constraint   = SpecificHeatDerivative(new_T_match,
                                                  &coef_high[0]);

  // generate the data to which the new lower branch is fit
  T_orig  = new double[num_points];
  Cp_orig = new double[num_points];

  for(int j=0; j<num_points; ++j) {
    T_orig[j]  = T_min + (new_T_match-T_min)*
      static_cast<double>(j)/(num_points-1.0);
    // set original Cp data according to the original T_match not the
    // new_T_match 
    if(T_orig[j] <= T_match) {
      Cp_orig[j] = SpecificHeat(T_orig[j],&coef_low[0]);
    } else {
      Cp_orig[j] = SpecificHeat(T_orig[j],&coef_high[0]);
    } 
  }

  // evaluate the new low temperature fits
  refit_success = PolynomialFit(num_points, // number of points generated to fit
                               4,          // polynomial degree
                               2,          // number of fixed point constraints
                               1,          // number of slope constraints
                               T_orig,
                               Cp_orig,
                               T_point_constraint,
                               Cp_point_constraint,
                               &T_slope_constraint,
                               &dCp_slope_constraint,
                               &refit_coef_low[0]);
  if(!refit_success) {
    printf("WARNING: PolynomialFit(...) failed,\n");
    printf("         retaining low temperature coefficients for Cp/R\n");
    for(int j=0; j<5; ++j) {
      refit_coef_low[j] = coef_low[j];
    }
  }
  // copy the original high coefficients
  for(int j=0; j<7; ++j) {
    refit_coef_high[j] = coef_high[j];
  }
  // copy the low coefficients for enthalpy and entropy
  refit_coef_low[5] = coef_low[5];
  refit_coef_low[6] = coef_low[6];

  // calculate the fixed temperature enthalpy and entropy from the original
  // thermo definition
  if(T_fixed < T_match) {
    H_fixed = Enthalpy(T_fixed,&coef_low[0]);
    S_fixed = Entropy( T_fixed,&coef_low[0]);
  } else {
    H_fixed = Enthalpy(T_fixed,&coef_high[0]);
    S_fixed = Entropy( T_fixed,&coef_high[0]);
  }
  // adjust the new thermo definition
  if(T_fixed < new_T_match) {
    // 1. adjust the low temperature coefficient branch
    MatchEnthalpy(T_fixed,H_fixed,&refit_coef_low[0]);
    MatchEntropy( T_fixed,S_fixed,&refit_coef_low[0]);
    // 2. adjust the high temprature coefficient branch to match
    H_fixed = Enthalpy(new_T_match,&refit_coef_low[0]);
    S_fixed = Entropy( new_T_match,&refit_coef_low[0]);

    MatchEnthalpy(new_T_match,H_fixed,&refit_coef_high[0]);
    MatchEntropy( new_T_match,S_fixed,&refit_coef_high[0]);
  } else {
    // 1. adjust the high temperature coefficient branch
    MatchEnthalpy(T_fixed,H_fixed,&refit_coef_high[0]);
    MatchEntropy( T_fixed,S_fixed,&refit_coef_high[0]);
    // 2. adjust the low temprature coefficient branch to match
    H_fixed = Enthalpy(new_T_match,&refit_coef_high[0]);
    S_fixed = Entropy( new_T_match,&refit_coef_high[0]);

    MatchEnthalpy(new_T_match,H_fixed,&refit_coef_low[0]);
    MatchEntropy( new_T_match,S_fixed,&refit_coef_low[0]);
  }

  delete [] T_orig;
  delete [] Cp_orig;
}


void RefitKeepHighGlobalTMatch(const int num_points,
                               const int T_resolution,
                               const double T_fixed,
                               const JanafThermoData *original,
                               const double global_T_match,
     		               JanafThermoData *refit)
{
  // copy the original min/max temperatures to the new JANAF fit, but
  // use the global T_match for the refit 
  refit->T_min   = original->T_min;
  refit->T_match = RoundToTenPower(-T_resolution, global_T_match);
  refit->T_max   = original->T_max;
  // initialize coefficients to zero
  for(int j=0; j<7; ++j) {
    refit->low_coef[j] = 0.0;
    refit->high_coef[j] = 0.0;
  }
  // check if only one temperature range is given
  if(refit->T_match <= refit->T_min) {
    // only in the high range, copy high range coefficients to both ranges
    for(int j=0; j<7; ++j) {
      refit->low_coef[j]  = original->high_coef[j];
      refit->high_coef[j] = original->high_coef[j];
    }
    return;
  } else if(refit->T_match >= refit->T_max) {
    // only in the low range, copy low range coefficients to both ranges
    for(int j=0; j<7; ++j) {
      refit->low_coef[j]  = original->low_coef[j];
      refit->high_coef[j] = original->low_coef[j];
    }
    return;
  }    
 
  RefitKeepHighGlobalTMatch(num_points,
                T_resolution, // decimal places for the match temperature
                T_fixed,
                original->T_min,
                original->T_match,
                original->T_max,
                original->low_coef,
                original->high_coef,
                refit->T_match,
                refit->low_coef,
                refit->high_coef);
}


void MatchEnthalpy(const double T_fixed,
                   const double H_fixed,
                   double coef[])

{
  double H_current = Enthalpy(T_fixed,&coef[0]);
  // a[5]/T coefficient
  coef[5] += (H_fixed - H_current)*T_fixed; 
}
void MatchEntropy(const double T_fixed,
                  const double S_fixed,
                  double coef[])

{
  double S_current = Entropy(T_fixed,&coef[0]);
  // a[6] coefficient
  coef[6] += S_fixed - S_current; 
}

double MinSpecificHeatDelta(const double T_tol,
                            const double T_min,
                            const double T_match,
                            const double T_max,
                            const double coef_low[],
                            const double coef_high[])
{
  const int MAX_ITER = 100;
  const double min_slope = 1.0e-32; // TODO using more meaniful bounds
                                    // e.g. the next iteration is out of range
                                    
  double T_current = T_match;
  double T_next;

  int iter = 0;
 
  if(fabs(SpecificHeatDerivativeDelta(T_current,coef_low,coef_high)) >
     min_slope) {
   T_next = T_current - 
    SpecificHeatDelta(T_current,coef_low,coef_high)/
    SpecificHeatDerivativeDelta(T_current,coef_low,coef_high);
  } else { // change in dCp/dT slope may cause nan
    printf("INFO: Near zero difference in dCp/dT detected at T = %.18g\n",
	   T_current);
    if(fabs(SpecificHeatDelta(T_current,coef_low,coef_high) > 
	    min_slope*(T_max-T_min))) {
       printf("WARNING: In MinSpecificHeatDelta(...),\n");
       printf("         failed to converge after %d iterations\n",iter);
       printf("         at T_match = %.18g\n", T_current);
       fflush(stdout);
       return 2.0*T_max;

    } else {
      return T_current;
    }
  }

  while(fabs(T_next-T_current) > T_tol && 
        T_min <= T_next &&
        T_next <= T_max) {

    ++iter;
    //printf("T_next = %12.8f\n",T_next);
    T_current = T_next;

    if(iter > MAX_ITER) {
      printf("WARNING: In MinSpecificHeatDelta(...),\n");
      printf("         failed to converge after %d iterations\n",iter);
      printf("         at T_match = %.18g\n", T_current);
      fflush(stdout);
      break;
    }
    if(fabs(SpecificHeatDerivativeDelta(T_current,coef_low,coef_high)) >
       min_slope) {
     T_next = T_current - 
      SpecificHeatDelta(T_current,coef_low,coef_high)/
      SpecificHeatDerivativeDelta(T_current,coef_low,coef_high);
    } else { // change in dCp/dT slope may cause nan
      printf("INFO: Near zero difference in dCp/dT detected at T = %.18g\n",
             T_current);
      if(fabs(SpecificHeatDelta(T_current,coef_low,coef_high) > 
         min_slope*(T_max-T_min))) {
        printf("WARNING: In MinSpecificHeatDelta(...),\n");
        printf("         failed to converge after %d iterations\n",iter);
        printf("         at T_match = %.18g\n", T_current);
        fflush(stdout);
        return 2.0*T_max;
      } else {
       return T_current;
      }
    }
    // Can lead to divide by zero
    // T_next = T_current - 
    //          SpecificHeatDelta(T_current,coef_low,coef_high)/
    //          SpecificHeatDerivativeDelta(T_current,coef_low,coef_high);
  } 
  if(iter > MAX_ITER) {
    return 2.0*T_max;
  }

  return T_next;
}

bool PolynomialFit(const int num_points,
                   const int degree,
                   const int num_fixed_points,
                   const int num_fixed_slopes,
                   const double x[],
                   const double f[],
                   const double x_fixed_points[],
                   const double f_fixed_points[],
                   const double x_fixed_slopes[],
                   const double dfdx_fixed_slopes[],
                   double coef[])
{
  // TODO: add dimension checking
  //  :
  //  :
  const int num_constraints = num_fixed_points+num_fixed_slopes;
  int num_rows_A=num_points;
  int num_cols_A=degree+1;
  int num_rows_B=num_constraints;
  int info;
  double *A, *B, *rhs, *constraint_rhs;
  int index;

  A   = new double[num_points*(degree+1)];
  rhs = new double[num_points];  

  // set up A matrix where A[k][j] = x[k]**j
  for(int k=0; k<num_points; ++k) {
    A[k] = 1.0;
  }
  index = num_points;
  for(int j=0; j<degree; ++j) {
    for(int k=0; k<num_points; ++k) {
      A[index]=A[index-num_points]*x[k];
      ++index;
    }
  }
  // copy the right-hand side to a temporary array
  for(int k=0; k<num_points; ++k) {
    rhs[k]=f[k];
  }

  if(num_fixed_points == 0 && num_fixed_slopes == 0) {
    // perform an unconstrained least squares optimization to
    // Solve for min |Ax-rhs|_2
    //info = LAPACKE_dgels(LAPACK_COL_MAJOR, // matrix element ordering
    //                     'N',              // Is matrix transposed?
    //                     num_rows_A,
    //                     num_cols_A,
    //                     1,                // number of right hand sides
    //                     A,                // matrix
    //                     num_rows_A,       // leading dimension of A matrix
    //                     rhs,              // right hand side matrix
    //                     num_rows_A);      // leading dimension of the
                                           // right hand side matrix
    info = LapackDGELS('N',             // Is matrix transposed?
                       num_rows_A,
                       num_cols_A,
                       1,                // number of right hand sides
                       A,                // matrix
                       num_rows_A,       // leading dimension of A matrix
                       rhs,              // right hand side matrix
                       num_rows_A);      // leading dimension of the
                                         // right hand side matrix

    if(info != 0) {
      //printf("WARNING: LAPACKE_dgels(...), returned info code = %d\n",
      //       info);
      printf("WARNING: LapackDGELS(...), returned info code = %d\n",
             info);
      delete [] A;
      delete [] rhs;
      return false;
    }
    // copy the solution to the coefficient array
    for(int j=0; j<=degree; ++j) {
      coef[j]=rhs[j];
    }

  } else {
    // allocate the constraint matrix and rhs
    B = new double[num_constraints*(degree+1)];
    constraint_rhs = new double[num_constraints];

    // build the constraint matrix for the fixed point constraints
    for(int j=0; j<=degree; ++j) {
      for(int k=0; k<num_fixed_points; ++k) {
        if(j==0) {
          // first column
          B[k+j*num_constraints] = 1.0;
        } else {
          B[k+j*num_constraints] =
            B[k+(j-1)*num_constraints]*x_fixed_points[k];
        }
      }
    }
    // build the constraint matrix for the fixed slopw constraints
    for(int j=0; j<=degree; ++j) {
      for(int k=num_fixed_points; k<num_constraints; ++k) {
        if(j==0) {
          // column 0
          B[k+j*num_constraints] = 0.0;
        } else if(j==1) {
          // column 1 
          B[k+j*num_constraints] = 1.0;
        }
        else {
          B[k+j*num_constraints] =
            static_cast<double>(j)*pow(x_fixed_slopes[k-num_fixed_points],
                                       j-1);
        }
      }
    }
    // build the constraint rhs
    for(int k=0; k<num_constraints; ++k) {
      if(k<num_fixed_points) {
        constraint_rhs[k] = f_fixed_points[k];
      } else {
        constraint_rhs[k] = dfdx_fixed_slopes[k-num_fixed_points]; 
      }
    }
    //info = LAPACKE_dgglse(LAPACK_COL_MAJOR, // matrix element ordering
    //                      num_rows_A,
    //                      num_cols_A,
    //                      num_rows_B,      // number of constraint rows
    //                      A,               // matrix
    //                      num_rows_A,      // leading dimension of A matrix
    //                      B,               // constraint matrix
    //                      num_rows_B,      // leading dimension of B matrix
    //                      rhs,             // right hand side vector
    //                      constraint_rhs,
    //                      coef);           // solution of linear constrained 
    //                                       // least squares
    info = LapackDGGLSE(num_rows_A,
                        num_cols_A,
                        num_rows_B,      // number of constraint rows
                        A,               // matrix
                        num_rows_A,      // leading dimension of A matrix
                        B,               // constraint matrix
                        num_rows_B,      // leading dimension of B matrix
                        rhs,             // right hand side vector
                        constraint_rhs,
                        coef);           // solution of linear constrained 
                                         // least squares
    if(info != 0) {
      //printf("WARNING: LAPACKE_dgglse(...), returned info code = %d\n",
      //       info);
      printf("WARNING: LapackDGGLSE(...), returned info code = %d\n",
             info);
      delete [] A;
      delete [] rhs;
      delete [] B;
      delete [] constraint_rhs;
      return false;
    }
    
    delete [] B;
    delete [] constraint_rhs;

  } // end constraint conditional

  delete [] A;
  delete [] rhs;
  return true;
}

double RoundToTenPower(const int n, const double a)
{
  double ten_power = pow(10.0,n);
  return ten_power*round(a/ten_power);
}
