
#include <cmath> //sqrt
#include <algorithm> //max

#include "utility_funcs.h"
#include "reactor_nvector_serial.h"

#include "nvector/nvector_serial.h"

ReactorNVectorSerial::ReactorNVectorSerial(std::shared_ptr<zerork::mechanism> mech_ptr) : 
  ReactorBase(),
  mech_ptr_(mech_ptr)
{
  num_species_ = mech_ptr_->getNumSpecies();
  num_variables_ = num_species_ + 1;
  num_steps_ = mech_ptr_->getNumSteps();

  sqrt_unit_round_ = sqrt(UNIT_ROUNDOFF);

  state_ = N_VNew_Serial(num_variables_);
  batch_mask_ = N_VNew_Serial(num_variables_);
  double* batch_mask_ptr = NV_DATA_S(batch_mask_);
  for(int k = 0; k < num_variables_; ++k) {
    batch_mask_ptr[k] = 1.0;
  }

  mol_wt_.resize(num_species_);
  inv_mol_wt_.resize(num_species_);
  net_production_rates_.resize(num_species_);
  energy_.resize(num_species_);
  cx_mass_.resize(num_species_);
  forward_rates_of_production_.resize(num_steps_);
  creation_rates_.resize(num_species_);
  destruction_rates_.resize(num_species_);
  concentrations_.resize(num_species_);

  jacobian_column_sums_.resize(num_variables_+1);
  jacobian_data_.clear();
  jacobian_row_indexes_.clear();

  preconditioner_column_sums_.resize(num_variables_+1);
  preconditioner_data_.clear();
  preconditioner_row_indexes_.clear();

  // set constant parameters
  mech_ptr_->getMolWtSpc(&mol_wt_[0]);
  for(int j=0; j < num_species_; ++j) {
    inv_mol_wt_[j]=1.0/mol_wt_[j];
  }

  SetupSparseJacobianArrays();

  weights_.assign(1,1.0);
}

ReactorNVectorSerial::~ReactorNVectorSerial()
{
  N_VDestroy(state_);
  N_VDestroy(batch_mask_);
}

N_Vector& ReactorNVectorSerial::GetStateNVectorRef() {
  return state_;
}

void ReactorNVectorSerial::SetBatchMaskNVector(int reactor_idx, N_Vector batch_mask) {
  assert(reactor_idx == 0);
  double* batch_mask_self_ptr = NV_DATA_S(batch_mask_);
  double* batch_mask_ptr = NV_DATA_S(batch_mask);
  memcpy(batch_mask_ptr, batch_mask_self_ptr, sizeof(double)*num_variables_);
}

void ReactorNVectorSerial::GetAbsoluteToleranceCorrection(N_Vector correction) {
  if(int_options_["abstol_dens"]) {
    double reactor_density = 1.0/inverse_density_;
    for(int j=0; j < num_species_; ++j) {
      double molar_density = reactor_density*inv_mol_wt_[j]*1.0e-3; //mks->cgs
      NV_Ith_S(correction,j) = 1.0/molar_density;
    }
    NV_Ith_S(correction,num_species_) = 1.0;
  } else {
    for(int j=0; j < num_variables_; ++j) {
      NV_Ith_S(correction,j) = 1.0;
    }
  }
}

void ReactorNVectorSerial::SetupSparseJacobianArrays() {
  const int dense_matrix_size = num_variables_*num_variables_;
  const int num_steps_=mech_ptr_->getNumSteps();
  std::vector<int> isNonZero(dense_matrix_size);

  jacobian_data_.assign(num_variables_+1,0);
  jacobian_diagonal_indexes_.assign(num_variables_,0);
  jacobian_last_row_indexes_.assign(num_variables_,0);

  num_noninteger_jacobian_nonzeros_ =
    mech_ptr_->getNonIntegerReactionNetwork()->GetNumJacobianNonzeros();

  noninteger_jacobian_.assign(num_noninteger_jacobian_nonzeros_,0.0);
  noninteger_sparse_id_.assign(num_noninteger_jacobian_nonzeros_,0);
  std::vector<int> noninteger_row_id;
  std::vector<int> noninteger_column_id;

  // initialize the dense nonzero flags
  for(int j=0; j<num_species_; ++j) {
      isNonZero[j*num_variables_+j]=1;    // mark the diagonal
      isNonZero[j*num_variables_+num_species_]=1; // mark the last row
  }
  for(int k=0; k<num_variables_; ++k) { // mark nSize rows in the last column
    isNonZero[num_species_*num_variables_+k]=1;
  }

  destruction_terms_.clear();
  creation_terms_.clear();
  // parse the system, filling in the Jacobian term data
  // Jacobian = d ydot(k)/ dy(j)
  for(int j=0; j<num_steps_; ++j) {
    int num_reactants=mech_ptr_->getOrderOfStep(j);
    int num_products=mech_ptr_->getNumProductsOfStep(j);
    for(int k=0; k<num_reactants; ++k) {
      int column_idx=mech_ptr_->getSpecIdxOfStepReactant(j,k); // species being perturbed
      // forward destruction
      for(int m=0; m < num_reactants; ++m) {
        int row_idx=mech_ptr_->getSpecIdxOfStepReactant(j,m); // species being destroyed
        isNonZero[column_idx*num_variables_+row_idx]=1; // mark location in densei

        reaction_indexes ri;
        ri.concentration_index = column_idx;
        ri.reaction_index = j;
        ri.sparse_index = row_idx;
        destruction_terms_.push_back(ri);
      }
      // forward creation
      for(int m=0; m < num_products; ++m) {
        int row_idx=mech_ptr_->getSpecIdxOfStepProduct(j,m); // species being created
        isNonZero[column_idx*num_variables_+row_idx]=1; // mark location in dense

        reaction_indexes ri;
        ri.concentration_index = column_idx;
        ri.reaction_index = j;
        ri.sparse_index = row_idx;
        creation_terms_.push_back(ri);
      }
    }
  }

  // non-integer reaction network
  if(num_noninteger_jacobian_nonzeros_ > 0) {

    noninteger_row_id.assign(num_noninteger_jacobian_nonzeros_, 0);
    noninteger_column_id.assign(num_noninteger_jacobian_nonzeros_, 0);

    mech_ptr_->getNonIntegerReactionNetwork()->GetJacobianPattern(
       &noninteger_row_id[0],&noninteger_column_id[0]);

    for(int j=0; j<num_noninteger_jacobian_nonzeros_; ++j) {

      int dense_id = noninteger_row_id[j]+noninteger_column_id[j]*num_variables_;
      isNonZero[dense_id]=1;

    }
  } // end if(num_noninteger_jacobian_nonzeros_ > 0)

  // count the number of nonzero terms in the dense matrix
  nnz_=1; // start the count at one so it can still serve as a flag
  jacobian_column_sums_[0]=0;
  for(int j=0; j<num_variables_; j++) {
    for(int k=0; k<num_variables_; k++) {
      if(isNonZero[j*num_variables_+k]==1) {
        isNonZero[j*num_variables_+k]=nnz_;
        nnz_ += 1;
      }
    }
    // after counting column j store the running total in column j+1
    // of the column sum
    jacobian_column_sums_[j+1]=nnz_-1;
  }
  // now at each nonzero term, isNonZero is storing the (address+1) in the
  // actual compressed column storage data array

  nnz_--; // decrement the count

  // allocate matrix data
  jacobian_data_.assign(nnz_,0.0);
  jacobian_row_indexes_.assign(nnz_,0);
  diagonal_indexes_.assign(num_variables_,0);
  last_row_indexes_.assign(num_variables_,0);

  // scan the the isNonZero array to determine the row indexes
  // and the special data addresses
  for(int j=0; j<num_variables_; ++j) {
    for(int k=0; k<num_variables_; ++k) {
      int nzAddr=isNonZero[j*num_variables_+k];
      if(nzAddr>0) {
        jacobian_row_indexes_[nzAddr-1]=k;
      }
    }
    // record the diagonal address
    diagonal_indexes_[j] = isNonZero[j*num_variables_+j]-1;
    last_row_indexes_[j] = isNonZero[(j+1)*num_variables_-1]-1;
  }

  // use the isNonZero array as a lookup to store the proper compressed
  // column data storage
  for(int j=0; j<destruction_terms_.size(); ++j) {
    int row_idx=destruction_terms_[j].sparse_index;
    int column_idx=destruction_terms_[j].concentration_index;
    int nzAddr=isNonZero[column_idx*num_variables_+row_idx];
    destruction_terms_[j].sparse_index=nzAddr-1; // reset to sparse addr
  }
  for(int j=0; j<creation_terms_.size(); ++j) {
    int row_idx=creation_terms_[j].sparse_index;
    int column_idx=creation_terms_[j].concentration_index;
    int nzAddr=isNonZero[column_idx*num_variables_+row_idx];
    creation_terms_[j].sparse_index=nzAddr-1; // reset to sparse addr
  }

  for(int j=0; j<num_noninteger_jacobian_nonzeros_; ++j) {
    int dense_id = noninteger_row_id[j]+noninteger_column_id[j]*num_variables_;
    int nzAddr=isNonZero[dense_id];
    noninteger_sparse_id_[j] = nzAddr-1;
  }
} 


#ifdef SUNDIALS2
int ReactorNVectorSerial::GetJacobianDense(long int N, double t, N_Vector y, N_Vector fy,
                                               DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  SetupJacobianSparse(t,y,fy,tmp1,tmp2,tmp3);
  int flag = SparseToDense(Jac);
  return flag;
}
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorNVectorSerial::GetJacobianDense(double t, N_Vector y, N_Vector fy,
                                               SUNMatrix Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  SetupJacobianSparse(t,y,fy,tmp1,tmp2,tmp3);
  int flag = SparseToDense(Jac);
  return flag;
}
#else
#error "Unsupported SUNDIALS version"
#endif


#ifdef SUNDIALS2
int ReactorNVectorSerial::SparseToDense(DlsMat Jac) {
  const int num_vars = GetNumStateVariables();
  for(int j=0; j < num_vars; ++j) {
    realtype* col = DENSE_COL(Jac,j);
    for(int k = jacobian_column_sums_[j]; k < jacobian_column_sums_[j+1]; ++k) {
      col[jacobian_row_indexes_[k]] = jacobian_data_[k];
    }
  }
  return 0;
}
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorNVectorSerial::SparseToDense(SUNMatrix Jac) {
  const int num_vars = GetNumStateVariables();
  for(int j=0; j < num_vars; ++j) {
    realtype* col = SUNDenseMatrix_Column(Jac,j);
    for(int k = jacobian_column_sums_[j]; k < jacobian_column_sums_[j+1]; ++k) {
      col[jacobian_row_indexes_[k]] = jacobian_data_[k];
    }
  }
  return 0;
}
#else
#error "Unsupported SUNDIALS version"
#endif


int ReactorNVectorSerial::SparseToDense(const std::vector<int> sums, const std::vector<int> idxs,
                                            const std::vector<double> vals, std::vector<double>* dense) {
  const int num_vars = GetNumStateVariables();
  dense->assign(num_vars*num_vars, 0.0);
  for(int j=0; j < num_vars; ++j) {
    for(int k = sums[j]; k < sums[j+1]; ++k) {
      (*dense)[j*num_vars + idxs[k]] = vals[k];
    }
  }
  return 0;
}

int ReactorNVectorSerial::JacobianSetup(double t, N_Vector y, N_Vector fy,
                                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
 if(int_options_["analytic"] == 1) {
   SetupJacobianSparse(t,y,fy,tmp1,tmp2,tmp3);
 } else {
   this->DividedDifferenceJacobian(t, y, fy, tmp1, tmp2, tmp3, &dense_jacobian_);
 }
 return 0;
}

int ReactorNVectorSerial::JacobianFactor(double gamma)
{
  bool dense_numeric = (int_options_["dense"] == 1) && (int_options_["analytic"] == 0);
  int mismatch = 0;
  if( !dense_numeric ) {
    mismatch = FormPreconditioner(gamma);
  }

  int flag = 0;
  if(int_options_["dense"] == 1) {
    if(int_options_["analytic"] == 1) {
      SparseToDense(preconditioner_column_sums_, preconditioner_row_indexes_, preconditioner_data_,
                    &dense_preconditioner_);
    } else {
      dense_preconditioner_.assign(num_variables_*num_variables_,0.0);
      for(int col = 0; col < num_variables_; ++col) {
        for(int row = 0; row < num_variables_; ++row) {
          const int idx = row*num_variables_ + col;
          dense_preconditioner_[idx] = -gamma*dense_jacobian_[idx];
        }
        dense_preconditioner_[col*num_variables_ + col] += 1.0;
      }
    }
    flag = lpm_.factor(num_variables_, num_variables_, dense_preconditioner_);
    dense_preconditioner_.clear();
  } else {
    //try refactor and fall back to full factor on fail
    //flag = 1;
    //if(mismatch == 0) {
      flag = slum_.refactor(preconditioner_data_);
    //}
    if(flag != 0) {
      flag = slum_.factor(preconditioner_row_indexes_,
                          preconditioner_column_sums_,
                          preconditioner_data_,
                          superlu_manager::CSC);
    }
  }
  return flag;
}

int ReactorNVectorSerial::FormPreconditioner(double gamma)
{
  int preconditioner_nnz = 0;
  double preconditioner_threshold = double_options_["preconditioner_threshold"];
  if(int_options_["iterative"] == 0) {
    preconditioner_threshold = 0.0;
  }

  std::vector<int> prev_preconditioner_column_sums = preconditioner_column_sums_;
  std::vector<int> prev_preconditioner_row_indexes = preconditioner_row_indexes_;

  preconditioner_data_.assign(nnz_,0.0);
  preconditioner_column_sums_.assign(num_variables_+1,0);
  preconditioner_row_indexes_.assign(nnz_,0);
  preconditioner_column_sums_[0] = 0;
  const int num_vars = num_variables_;

  if(int_options_["analytic"] == 0) {
    for(int j=0; j<num_vars; ++j) {
      preconditioner_column_sums_[j+1] = preconditioner_column_sums_[j];
      for(int k=jacobian_column_sums_[j]; k<jacobian_column_sums_[j+1]; ++k) {
        int row = jacobian_row_indexes_[k];
        bool diag = row == j;
        double element_value = -gamma*dense_jacobian_[j*num_variables_ + row];
        if(diag) {
          element_value += 1.0;
        }

        if(diag || fabs(element_value) > preconditioner_threshold) { // directly copy diagonal and any 
          //off-diagonal terms larger than  tol
          preconditioner_data_[preconditioner_nnz]=element_value;
          preconditioner_row_indexes_[preconditioner_nnz] = row;
          preconditioner_column_sums_[j+1] += 1;
          preconditioner_nnz += 1;
        }
      }
    }
  } else {
    for(int j=0; j<num_vars; ++j) {
      preconditioner_column_sums_[j+1] = preconditioner_column_sums_[j];
      for(int k=jacobian_column_sums_[j]; k<jacobian_column_sums_[j+1]; ++k) {
        int row = jacobian_row_indexes_[k];
        bool diag = row == j;
        double element_value = -gamma*jacobian_data_[k];
        if(diag) {
          element_value += 1.0;
        }

        if(diag || fabs(element_value) > preconditioner_threshold) { // directly copy diagonal and any 
          //off-diagonal terms larger than  tol
          preconditioner_data_[preconditioner_nnz]=element_value;
          preconditioner_row_indexes_[preconditioner_nnz] = row;
          preconditioner_column_sums_[j+1] += 1;
          preconditioner_nnz += 1;
        }
      }
    }
  }

  //Check for mismatch
  for(int j=1; j<=num_vars; ++j) {
    if(prev_preconditioner_column_sums[j] != preconditioner_column_sums_[j]) {
      return 1;
    }
  }
  for(int j=0; j<preconditioner_nnz; ++j) {
    if(prev_preconditioner_row_indexes[j] != preconditioner_row_indexes_[j]) {
      return 1;
    }
  }
  return 0;
}

int ReactorNVectorSerial::JacobianSolve(double t, N_Vector y, N_Vector fy,
                                            N_Vector r, N_Vector z, N_Vector tmp)
{
  double * r_ptr = NV_DATA_S(r);
  double * z_ptr = NV_DATA_S(z);
  int flag = 0;
  if(int_options_["dense"] == 1) {
    flag = lpm_.solve(num_variables_, r_ptr, z_ptr);
  } else {
    flag = slum_.solve(num_variables_, r_ptr, z_ptr);
  }
  return flag;
}


int ReactorNVectorSerial::SetupJacobianSparse(realtype t, N_Vector y,N_Vector fy,
                                                  N_Vector tmp1,N_Vector tmp2,N_Vector tmp3)
{
  double *y_ptr = NV_DATA_S(y);
  double *fy_ptr = NV_DATA_S(fy);
  double *tmp1_ptr = NV_DATA_S(tmp1);
  double *tmp2_ptr = NV_DATA_S(tmp2);
  double *tmp3_ptr = NV_DATA_S(tmp3);

  //double startTime = getHighResolutionTime();

  const int num_vars = num_variables_;
  const int num_spec = num_species_;
  const int num_destruction_terms = destruction_terms_.size();
  const int num_creation_terms = creation_terms_.size();
  const double RuTemp= mech_ptr_->getGasConstant() * y_ptr[num_spec];
  const double min_mass_frac = double_options_["min_mass_fraction"];

  // set tmp1 to the strictly positive mass fraction array
  // set tmp2 to 1/C where C is the concentration for the strictly positive
  // mass fraction array
  for(int j=0; j < num_spec; j++) {
    tmp1_ptr[j]=y_ptr[j];
    if(fabs(tmp1_ptr[j]) < min_mass_frac) {
      if(tmp1_ptr[j] < 0.0) {
        tmp1_ptr[j]= -min_mass_frac;
      } else {
        tmp1_ptr[j]=  min_mass_frac;
      }
    }
    tmp2_ptr[j]=mol_wt_[j]/tmp1_ptr[j]*inverse_density_;
  }

  tmp1_ptr[num_spec]=y_ptr[num_spec];

  // calculate the reaction info for the strictly positive case
  //double start_time_deriv = getHighResolutionTime();
  GetTimeDerivative(t,tmp1,tmp3);
  //double deriv_time = getHighResolutionTime() - start_time_deriv;

  // set the full sparse array
  jacobian_data_.assign(nnz_,0.0);

  // process the forward destruction terms
  for(int j=0; j < num_destruction_terms; ++j) {
      int conc_idx = destruction_terms_[j].concentration_index;
      int rxn_idx    = destruction_terms_[j].reaction_index;
      int sparse_idx = destruction_terms_[j].sparse_index;
      jacobian_data_[sparse_idx]-=forward_rates_of_production_[rxn_idx]*tmp2_ptr[conc_idx];
  }

  // process the forward creation terms
  for(int j=0; j < num_creation_terms; ++j) {
      int conc_idx = creation_terms_[j].concentration_index;
      int rxn_idx    = creation_terms_[j].reaction_index;
      int sparse_idx = creation_terms_[j].sparse_index;
      jacobian_data_[sparse_idx]+=forward_rates_of_production_[rxn_idx]*tmp2_ptr[conc_idx];
  }

  // process the non-integer Jacobian information
  const int num_noninteger_jacobian_nonzeros =
    num_noninteger_jacobian_nonzeros_;

  for(int j=0; j<num_noninteger_jacobian_nonzeros; ++j) {
    noninteger_jacobian_[j] = 0.0;
  }

  mech_ptr_->getNonIntegerReactionNetwork()->GetSpeciesJacobian(
    tmp2_ptr,
    &forward_rates_of_production_[0],
    &noninteger_jacobian_[0]);

  for(int j=0; j<num_noninteger_jacobian_nonzeros; ++j) {
    jacobian_data_[noninteger_sparse_id_[j]] += noninteger_jacobian_[j];
  }

  // At this point sMptr stores d(wdot[k])/dC[j] ignoring the contribution
  // of perturbations in the third body species

  // ---------------------------------------------------------------------
  // compute d(Tdot)/dy[j]

  // SPARSE:
  // step 1: compute matrix vector product (d(wdot[k])/dC[j])*E[k]
  for(int j=0; j<num_spec; j++) { // column number
    int last_row_idx = last_row_indexes_[j];

    jacobian_data_[last_row_idx]=0.0;
    //N.B. skip last row in this loop
    for(int k=jacobian_column_sums_[j]; k<jacobian_column_sums_[j+1]-1; k++) {
      int row=jacobian_row_indexes_[k];
      jacobian_data_[last_row_idx]+=jacobian_data_[k]*energy_[row];
    }
    jacobian_data_[last_row_idx] *= RuTemp;
  }

  // step 2: add the Tdot*Cv[j] term to d(Tdot)/dy[j]
  for(int j=0; j<num_spec; j++) {
    int last_row_idx = last_row_indexes_[j];
    //TODO: Tdot! (with non-zero concentrations) (CHECK THIS)
    jacobian_data_[last_row_idx]+=tmp3_ptr[num_spec]*cx_mass_[j]*mol_wt_[j];
  }

  // step 3: divide by -1/(molWt[j]*meanCvMass)
  double multFact = -1.0/(mean_cx_mass_);
  for(int j=0; j<num_spec; j++) {
    int last_row_idx = last_row_indexes_[j];
    jacobian_data_[last_row_idx] *= inv_mol_wt_[j]*multFact;
  }

  //Convert concentration to mass fraction
  for(int j=0; j<num_spec; j++) { // jth column
    for(int k=jacobian_column_sums_[j]; k<jacobian_column_sums_[j+1]-1; k++) {
      int row=jacobian_row_indexes_[k];
      jacobian_data_[k] *= mol_wt_[row]*inv_mol_wt_[j];
    }
  }

  // At this point Mtx stores d(Tdot[k])/dy[j] ignoring the contribution
  // of perturbations in the third body species

  // --------------------------------------------------------------------- 
  // calculate d(ydot[k])/dT

  // step 1: perturb the temperature
  double delta_temp=y_ptr[num_spec]*sqrt_unit_round_;
  tmp1_ptr[num_spec]+=delta_temp;
  delta_temp=tmp1_ptr[num_spec]-y_ptr[num_spec];
  double multiplier = 1.0/delta_temp;

  for(int j=0; j<num_spec; ++j) {
    tmp1_ptr[j]=y_ptr[j];
  }

  // step 2: calculate ydot at Temp+dTemp
  //start_time_deriv = getHighResolutionTime();
  GetTimeDerivative(t,tmp1,tmp3);
  //deriv_time += getHighResolutionTime() - start_time_deriv;

  // step 3: approximate d(ydot[k])/dT with finite difference
  for(int k=jacobian_column_sums_[num_spec]; k<jacobian_column_sums_[num_vars]; ++k) {
    int row=jacobian_row_indexes_[k];
    jacobian_data_[k]=(tmp3_ptr[row]-fy_ptr[row])*multiplier;
  }

  //jacobian_setup_time += getHighResolutionTime() - start_time - deriv_time;
  //jacobian_setups += 1;
  //bool print_jacobian = false;
  //if(print_jacobian) {
  //  print_sp_matrix(num_variables_, num_variables_, &jacobian_column_sums_[0], &jacobian_row_indexes_[0], &jacobian_data_[0]);
  //  std::vector<double> dense;
  //  SparseToDense(jacobian_column_sums_, jacobian_row_indexes_, jacobian_data_, &dense);
  //  printf("\n");
  //}
  return 0;
}


int ReactorNVectorSerial::RootFunction(double t, N_Vector y, double *root_function)
{
  double ignition_temperature = initial_temperature_ + double_options_["delta_temperature_ignition"];
  double current_temperature = NV_Ith_S(y,num_species_)*double_options_["reference_temperature"];
  root_function[0] = ignition_temperature - current_temperature;
  return 0;
}


int ReactorNVectorSerial::GetNumStateVariables()
{
  return num_variables_;
}


int ReactorNVectorSerial::GetNumRootFunctions()
{
  return 1;
}


void ReactorNVectorSerial::DividedDifferenceJacobian(double t, N_Vector y, N_Vector fy,
                                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3,
                                                     std::vector<double>* dense_jacobian) {
  static const int n = num_variables_;
  static const double uround = 1.0e-16;
  double* yd = NV_DATA_S(y);
  double* dy0d = NV_DATA_S(tmp1);
  double* dy1d = NV_DATA_S(tmp2);

  this->GetTimeDerivative(t,y,tmp1);
  dense_jacobian->assign(num_variables_*num_variables_,0.0);
  for(int i = 0; i < n; ++i) {
     double ysafe=yd[i];
     double delta=sqrt(uround*std::max(1.0e-5,fabs(ysafe)));
     yd[i]=ysafe+delta;
     this->GetTimeDerivative(t,y,tmp2);
     for(int j = 0; j < n; ++j) {
       (*dense_jacobian)[i*n + j] = (dy1d[j] - dy0d[j])/delta;
     }
     yd[i] = ysafe;
  }
}


