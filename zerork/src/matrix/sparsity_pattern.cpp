#include <stdlib.h>
#include <string.h>

#include <limits>

#include "sparsity_pattern.h"

namespace zerork {

SparsityPattern::SparsityPattern(mechanism *mechp, const char *type)
{
  char *jacobian_type;
  int type_length = strlen(type)+1;

  jacobian_type = new char[type_length];
  strcpy(jacobian_type,type);
  // convert jacobian_type to lower case
  for(int j=0; j<type_length; ++j) {
    if('A' <= jacobian_type[j] && jacobian_type[j] <= 'Z') {
      jacobian_type[j]+=('a' - 'A');
    }
  }

  if(strcmp(jacobian_type,"cmt")==0) {
    // Create Jacobian sparsity pattern for a state vector consisting of
    // species [C]oncentrations, total [M]ixture concentration, and
    // [T]emperature.
    num_equations_ = mechp->getNumSpecies();

    AddBasicStoichiometry(mechp);
    AddEnhancedSpecies(mechp);
    AddMixtureConcentration(mechp);
    AddFullRowColumn(); // for temperature
    AddDiagonal();

  } else {

    printf("ERROR: In SparsityPattern(zerork::mechanism *, const char *),\n");
    printf("       Jacobian type \"%s\" not recognized.\n",type);
    printf("       Current supported type(s): \"CMT\".\n");
    exit(-1);
  }

  CreatePositionAddress();

  delete [] jacobian_type;
}

// GetCompressedColumnStorage(int row_id[], int col_sum[])
// 
// Return the sparsity pattern as represented by compressed column storage
// to the two array pointers row_id (row index) and col_sum (column sum).
//
// length(row_id)  = number of nonzero terms
// length(col_sum) = number of equations plus one
//
// return value is the number of nonzero terms.
int SparsityPattern::GetCompressedColumnStorage(int row_id[],
                                                int col_sum[])
{
  size_t nnz = 0;
  
  col_sum[0] = 0;

  for(int j=0; j<num_equations_; ++j) {   // column j
    col_sum[j+1]=col_sum[j];
    for(int k=0; k<num_equations_; ++k) { // row k
      if(PositionCount(k,j) > 0) {
        row_id[nnz] = k;
	++col_sum[j+1];
        ++nnz;
      }
    }
  }
  if(nnz != position_count_.size()) {
    printf("ERROR: In SparsityPattern::GetCompressedColumnStorage(...),\n");
    printf("       number of nonzero element addresses recorded %lu\n",
           nnz);
    printf("       does not equal the size the position count map %lu.\n",
           position_count_.size());
    exit(-1);
  }
  return position_count_.size();
}

void SparsityPattern::AddBasicStoichiometry(mechanism *mechp)
{
  const int num_steps_const = mechp->getNumSteps(); 

  for(int j=0; j<num_steps_const; ++j) {

    int num_reactants = mechp->getOrderOfStep(j);
    int num_products  = mechp->getNumProductsOfStep(j);
    for(int k=0; k<num_reactants; ++k) {
      // derivative w.r.t. to the k^th reactant
      // affects the net rate of the m^th reactant
      for(int m=0; m<num_reactants; ++m) {
        // row (species rate)
        // col (species deriv)
        IncrementPosition(mechp->getSpecIdxOfStepReactant(j,m),
                          mechp->getSpecIdxOfStepReactant(j,k));
      }
      // affects the net rate of the m^th product
      for(int m=0; m<num_products; ++m) {
        // row (species rate)
        // col (species deriv)
        IncrementPosition(mechp->getSpecIdxOfStepProduct(j,m),
                          mechp->getSpecIdxOfStepReactant(j,k));
      }
    } // k-loop (num_reactants)
  } // j-loop (num_steps_const)
}


void SparsityPattern::AddDiagonal()
{
  const int num_equations_const = num_equations_;

  for(int j=0; j<num_equations_const; ++j) {
    IncrementPosition(j,j);
  }
}
void SparsityPattern::AddFullRowColumn()
{
  const int num_equations_const = num_equations_;
  
  for(int j=0; j<num_equations_const; ++j) {
    IncrementPosition(num_equations_const,j); // fill in a new row
  }  
  ++num_equations_;
  for(int j=0; j<num_equations_const+1; ++j) {
    IncrementPosition(j,num_equations_const); // fill in a new column
  }  

}

// Add new row and column to the Jacobian for the total mixture concentration.
void SparsityPattern::AddMixtureConcentration(mechanism *mechp)
{
  const int num_equations_const = num_equations_;
  const int num_reactions_const = mechp->getNumReactions();

  for(int j=0; j<num_equations_const; ++j) {
    IncrementPosition(num_equations_const,j); // fill in a new row
  }  
  ++num_equations_;
  // new column contains non-zero elements for any reactant involved in 
  // a third body or falloff reaction
  for(int j=0; j<num_reactions_const; ++j) {
    if(mechp->isThirdBodyReaction(j) == 1 || 
       mechp->isThirdBodyFalloffReaction(j) == 1) {

      int fwd_step_id = mechp->getStepIdxOfRxn(j, 1);
      for(int k=0; k<mechp->getOrderOfStep(fwd_step_id); ++k) {
        IncrementPosition(mechp->getSpecIdxOfStepReactant(fwd_step_id,k),
                          num_equations_const);
      }
      for(int k=0; k<mechp->getNumProductsOfStep(fwd_step_id); ++k) {
        IncrementPosition(mechp->getSpecIdxOfStepProduct(fwd_step_id,k),
                          num_equations_const);
      }
    }
  }
}

void SparsityPattern::AddEnhancedSpecies(mechanism *mechp)
{
  const int num_reactions_const = mechp->getNumReactions();
  std::vector<int> enhanced_species_id;
  std::vector<double> enhanced_species_alpha;

  for(int j=0; j<num_reactions_const; ++j) {
 
    if(mechp->getNumEnhancedSpeciesOfStep(j) > 0) {
      
      int fwd_step_id = mechp->getStepIdxOfRxn(j, 1);

      mechp->getEnhancementFactorsOfStep(j,
                                         &enhanced_species_id,
                                         &enhanced_species_alpha);
      // Add the derivative of each enhanced third body species for the net
      // rate of each reactant and product.
      for(int k=0; k<mechp->getNumEnhancedSpeciesOfStep(j); ++k) {
        
        for(int m=0; m<mechp->getOrderOfStep(fwd_step_id); ++m) {
          IncrementPosition(mechp->getSpecIdxOfStepReactant(fwd_step_id,m),
                            enhanced_species_id[k]);
        }
        for(int m=0; m<mechp->getNumProductsOfStep(fwd_step_id); ++m) {
          IncrementPosition(mechp->getSpecIdxOfStepProduct(fwd_step_id,m),
                            enhanced_species_id[k]);
        }
      }
    } // end if(mechp->getNumEnhancedSpeciesOfStep(j) > 0)
  } // end for(int j=0; j<num_reactions_const; ++j)
}
  
int SparsityPattern::PositionCount(const int row_id, const int col_id)
{
  std::pair<int,int> search_position;
  std::map<std::pair<int,int>, int>::iterator iter;

  search_position = std::make_pair(row_id,col_id);
  iter = position_count_.find(search_position);

  if(iter != position_count_.end()) {
    return iter->second;
  } else {
    return 0;
  }
}
 
int SparsityPattern::PositionAddress(const int row_id, const int col_id)
{
  std::pair<int,int> search_position;
  std::map<std::pair<int,int>, int>::iterator iter;

  search_position = std::make_pair(row_id,col_id);
  iter = position_address_.find(search_position);

  if(iter != position_address_.end()) {
    return iter->second;
  } else {
    return std::numeric_limits<int>::min();
  }
}

void  SparsityPattern::IncrementPosition(const int row_id, const int col_id)
{
  std::pair<int,int> search_position;
  std::pair<std::pair<int,int>,int> new_position;

  std::map<std::pair<int,int>, int>::iterator iter;

  search_position = std::make_pair(row_id,col_id);
  iter = position_count_.find(search_position);

  if(iter != position_count_.end()) {
    iter->second++;
  } else {
    new_position.first  = search_position; // map-key
    new_position.second = 1;               // map value
    position_count_.insert(new_position); 
  }
}

void SparsityPattern::CreatePositionAddress()
{
  int nnz=0;
  std::pair<int,int> search_position;
  std::pair<std::pair<int,int>,int> new_address;

  for(int j=0; j<num_equations_; ++j) {   // column j
    for(int k=0; k<num_equations_; ++k) { // row k
      if(PositionCount(k,j) > 0) {
        search_position = std::make_pair(k,j);
        new_address.first  = search_position;  // map-key
        new_address.second = nnz;
        position_address_.insert(new_address);
        // note map::insert does not insert duplicates
        ++nnz;
      }
    }
  }
  if(nnz != static_cast<int>(position_address_.size())) {
    printf("ERROR: In SparsityPattern::CreatePositionAddress(...),\n");
    printf("       number of nonzero element addresses recorded %d\n",
           nnz);
    printf("       does not equal the size the position address map %lu.\n",
           position_address_.size());
    exit(-1);
  }
}

}// namespace zerork
