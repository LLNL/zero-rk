#ifndef ZERORK_REACTOR_MANAGER_H
#define ZERORK_REACTOR_MANAGER_H


#include <string>
#include <vector>
#include <memory>
#include <fstream>

#include "zerork_reactor_manager_base.h"

#include "reactor_base.h"

#include "zerork/mechanism.h"

class ZeroRKReactorManager : public ZeroRKReactorManagerBase
{
 public:
  ZeroRKReactorManager(const char *input_filename,
                       const char* mech_filename,
                       const char* therm_filename);
  virtual ~ZeroRKReactorManager() {};

  void SetInputVariables(int n_cycle,
                         double time,
                         double dt,
                         int n_reactors, 
                         double* T,
                         double* P,
                         double* mass_fractions);

  void SetAuxFieldPointer(zerork_field_type ft, double* field_pointer);

  void LoadBalance();
  void SolveReactors();
  void RedistributeResults();
  void PostSolve();

 private:
  std::shared_ptr<zerork::mechanism> mech_ptr_;

  int n_reactors_self_;
  int n_reactors_other_;

  int num_species_;
  int num_species_stride_;
  int rank_;
  int sorted_rank_;
  int nranks_;
  int root_rank_;
  int load_balance_;

  int tx_count_per_reactor_;
  std::vector<int> comm_mtx_row_sum_;
  std::vector<int> comm_mtx_col_idx_;
  std::vector<int> comm_mtx_reactors_;
  std::vector<size_t> sorted_deficit_idxs_;
  std::vector<size_t> sorted_reactor_idxs_;

  double dt_calc_;
  double* T_self_;
  double* P_self_;
  double* mf_self_;
  double* dpdt_self_;
  double* rc_self_;
  double* rg_self_;
  double* y_src_self_;
  double* e_src_self_;
  int* nstep_self_;

  std::vector<double> T_other_;
  std::vector<double> P_other_;
  std::vector<double> mf_other_;
  std::vector<double> dpdt_other_;
  std::vector<double> e_src_other_;
  std::vector<double> y_src_other_;
  std::vector<double> rc_other_;
  std::vector<double> rg_other_;
  std::vector<int> nstep_other_;

  bool dpdt_owned_;
  bool rc_owned_;
  bool rg_owned_;
  bool y_src_defined_;
  bool e_src_defined_;
  std::vector<double> dpdt_default_;
  std::vector<double> rc_default_;
  std::vector<double> rg_default_;

  int n_calls_;
  int n_cycle_;
  int n_solve_;
  int n_steps_;
  int n_weight_updates_;
  int n_reactors_min_;
  int n_reactors_max_;
  double sum_reactor_time_;
  std::vector<double> rank_weights_;
  std::vector<int> n_reactors_solved_ranks_;
  std::vector<double> all_time_ranks_;
  std::ofstream reactor_log_file_;

  std::unique_ptr<ReactorBase> reactor_ptr_;

  static const int EXCHANGE_SEND_TAG_ = 42;
  static const int EXCHANGE_RETURN_TAG_ = 43;

  int RecvReactors(size_t send_rank);
  void SendReactors(std::vector<size_t> send_reactor_idxs, size_t recv_rank);

  void UpdateRankWeights();
  void ProcessPerformance();
};

#endif
