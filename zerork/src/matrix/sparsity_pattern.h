#ifndef SPARSITY_PATTERN_H
#define SPARSITY_PATTERN_H

#include "zerork/mechanism.h"

namespace zerork {

class SparsityPattern
{
 public:
  SparsityPattern(mechanism *mechp, const char *type);
  int GetCompressedColumnStorage(int row_id[], int col_sum[]);
  int PositionCount(const int row_id, const int col_id);
  int PositionAddress(const int row_id, const int col_id);
  int num_non_zeros() const {return position_count_.size();}
  int num_equations() const {return num_equations_;}  

 private:
  int num_equations_;
  std::map<std::pair<int,int>, int> position_count_;
  std::map<std::pair<int,int>, int> position_address_;

  void IncrementPosition(const int row_id, const int col_id);
  void AddFullRowColumn();
  void AddBasicStoichiometry(mechanism *mechp);
  void AddEnhancedSpecies(mechanism *mechp);
  void AddMixtureConcentration(mechanism *mechp);
  void AddDiagonal();
  void CreatePositionAddress();
};

}// namespace zerork

#endif

