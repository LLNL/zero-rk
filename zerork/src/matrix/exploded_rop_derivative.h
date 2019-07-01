#ifndef EXPLODED_ROP_DERIVATIVE_H_
#define EXPLODED_ROP_DERIVATIVE_H_

#include "zerork/mechanism.h"
#include "sparsity_pattern.h"

namespace zerork {

class ExplodedROPDerivative 
{
 public:
  ExplodedROPDerivative(mechanism *mechp,
                        SparsityPattern *pattern);
  ~ExplodedROPDerivative();

  int num_creation_terms() const {return num_creation_terms_;}
  int num_destruction_terms() const {return num_destruction_terms_;}
  void GetSpeciesJacobian(const double inv_conc[], 
                          const double step_rop[],
                          double jacobian[]) const;

 private:
  typedef struct{
    int row_id;
    int col_id;
    int step_id;
    int sparse_id; 
  } SparseJacobianLocator;

  int num_creation_terms_;
  int num_destruction_terms_;
  int num_non_zeros_;

  std::vector<SparseJacobianLocator> creation_list_;
  std::vector<SparseJacobianLocator> destruction_list_;

  void BuildROPLists(mechanism *mechp,
                     SparsityPattern *pattern);  
};

} // namespace zerork
#endif
