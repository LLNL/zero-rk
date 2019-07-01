#include "exploded_rop_derivative.h"

namespace zerork {

ExplodedROPDerivative::ExplodedROPDerivative(mechanism *mechp,
                                             SparsityPattern *pattern)
{
  BuildROPLists(mechp,pattern);
  num_non_zeros_ = pattern->num_non_zeros();
  // Could add a sorting step to the of the lists to possibly improve 
  // cache access.

  // Could add a copy operation to replace vectors with arrays.
}

ExplodedROPDerivative::~ExplodedROPDerivative()
{0;}

void ExplodedROPDerivative::BuildROPLists(mechanism *mechp,
                                          SparsityPattern *pattern)
{
  const int num_steps = mechp->getNumSteps();
  int num_reactants,num_products;
  int row_id,col_id;
  SparseJacobianLocator new_term;

  for(int j=0; j<num_steps; ++j) {
    num_reactants = mechp->getOrderOfStep(j);
    num_products  = mechp->getNumProductsOfStep(j);

    for(int k=0; k<num_reactants; ++k) {
      // col_id = species being perturbed or the species of which the partial
      // derivative is to be taken of the step's rate of progress
      col_id = mechp->getSpecIdxOfStepReactant(j,k);
      
      // species destruction terms (i.e. other reactants in the reaction step)
      for(int m=0; m<num_reactants; ++m) {
        // row_id = species whose destruction rate is impacted by a
        // perturbation in species col_id's concentration
        row_id = mechp->getSpecIdxOfStepReactant(j,m);

        new_term.row_id    = row_id;
        new_term.col_id    = col_id;
        new_term.step_id   = j;
        new_term.sparse_id = pattern->PositionAddress(row_id,col_id);

        destruction_list_.push_back(new_term);
      }
      // species creation terms (i.e. products in the reaction step)
      for(int m=0; m<num_products; ++m) {
        // row_id = species whose destruction rate is impacted by a
        // perturbation in species col_id's concentration
        row_id = mechp->getSpecIdxOfStepProduct(j,m);

        new_term.row_id    = row_id;
        new_term.col_id    = col_id;
        new_term.step_id   = j;
        new_term.sparse_id = pattern->PositionAddress(row_id,col_id);

        creation_list_.push_back(new_term);
      }
    } // end for(int k=0; k<num_reactants; ++k)
  } // end for(int j=0; j<num_steps; ++j)

  num_creation_terms_    = creation_list_.size();
  num_destruction_terms_ = destruction_list_.size();
}

void ExplodedROPDerivative::GetSpeciesJacobian(const double inv_conc[], 
                                               const double step_rop[],
                                               double jacobian[]) const
{
  const int num_non_zeros_const         = num_non_zeros_;
  const int num_creation_terms_const    = num_creation_terms_;
  const int num_destruction_terms_const = num_destruction_terms_;

  for(int j=0; j<num_non_zeros_const; ++j) {
    jacobian[j]=0.0;
  }
  for(int j=0; j<num_creation_terms_const; ++j) {
    jacobian[creation_list_[j].sparse_id] += 
      (step_rop[creation_list_[j].step_id]*
       inv_conc[creation_list_[j].col_id]);
  }
  for(int j=0; j<num_destruction_terms_const; ++j) {
    jacobian[destruction_list_[j].sparse_id] -= 
      (step_rop[destruction_list_[j].step_id]*
       inv_conc[destruction_list_[j].col_id]);
  }
}


} // namespace zerork

