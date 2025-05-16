#ifndef ZERORK_REACTOR_MANAGER_H
#define ZERORK_REACTOR_MANAGER_H


#include <string>
#include <vector>
#include <memory>
#include <fstream>

#include "zerork_reactor_manager_base.h"

#include "reactor_base.h"

#include "zerork/mechanism.h"
#ifdef ZERORK_GPU
#include "zerork/mechanism_cuda.h"
#endif

class ZeroRKReactorManager : public ZeroRKReactorManagerBase
{
 public:
  ZeroRKReactorManager();
  virtual ~ZeroRKReactorManager() {};

  zerork_status_t ReadOptionsFile(const std::string& options_filename);
  zerork_status_t LoadMechanism();

  zerork_status_t SetInputVariables(int n_cycle,
                         double time,
                         double dt,
                         int n_reactors,
                         double* T,
                         double* P,
                         double* mass_fractions);

  zerork_status_t SetAuxFieldPointer(zerork_field_t ft, double* field_pointer);
  zerork_status_t SetCallbackFunction(zerork_callback_fn fn, void* cb_fn_data);
  zerork_status_t SetReactorIDs(int* reactor_ids);

  zerork_status_t FinishInit();
  zerork_status_t LoadBalance();
  zerork_status_t SolveReactors();
  zerork_status_t RedistributeResults();
  zerork_status_t PostSolve();

 private:
  std::shared_ptr<zerork::mechanism> mech_ptr_;
#ifdef ZERORK_GPU
  std::shared_ptr<zerork::mechanism_cuda> mech_cuda_ptr_;
#endif

  int n_reactors_self_;
  int n_reactors_other_;

  int num_species_;
  int num_species_stride_;
  int rank_;
  int sorted_rank_;
  int nranks_;
  int root_rank_;

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
  double* root_times_self_;
  double* temp_delta_self_;
  int* nstep_self_;
  int* reactor_ids_self_;

  std::vector<double> T_other_;
  std::vector<double> P_other_;
  std::vector<double> mf_other_;
  std::vector<double> dpdt_other_;
  std::vector<double> e_src_other_;
  std::vector<double> y_src_other_;
  std::vector<double> rc_other_;
  std::vector<double> rg_other_;
  std::vector<double> root_times_other_;
  std::vector<double> temp_delta_other_;
  std::vector<int> nstep_other_;
  std::vector<int> reactor_ids_other_;

  bool rc_owned_;
  bool rg_owned_;
  bool root_times_owned_;
  bool temp_delta_owned_;
  bool dpdt_defined_;
  bool e_src_defined_;
  bool y_src_defined_;
  bool reactor_ids_defined_;
  std::vector<double> rc_default_;
  std::vector<double> rg_default_;
  std::vector<double> root_times_default_;
  std::vector<double> temp_delta_default_;

  bool tried_init_;
  int n_calls_;
  int n_cycle_;
  int n_gpu_solve_;
  int n_gpu_solve_no_temperature_;
  int n_cpu_solve_;
  int n_cpu_solve_no_temperature_;
  int n_steps_cpu_;
  int n_steps_gpu_;
  int n_gpu_groups_;
  int n_weight_updates_;
  double gpu_multiplier_;
  double sum_cpu_reactor_time_;
  double sum_gpu_reactor_time_;
  double avg_reactor_time_;
  std::vector<double> rank_weights_;
  std::vector<int> n_reactors_solved_ranks_;
  std::vector<double> all_time_ranks_;
  std::ofstream reactor_log_file_;
  zerork_callback_fn cb_fn_;
  void* cb_fn_data_;

  std::unique_ptr<ReactorBase> reactor_ptr_;
#ifdef ZERORK_GPU
  int n_cpu_ranks_;
  int n_gpu_ranks_;
  int gpu_id_;
  std::vector<int> rank_has_gpu_;
  std::unique_ptr<ReactorBase> reactor_gpu_ptr_;
  void AssignGpuId();
  void UpdateRankWeights();
#endif

  static const int EXCHANGE_SEND_TAG_ = 42;
  static const int EXCHANGE_RETURN_TAG_ = 43;

#if USE_MPI
  int RecvReactors(size_t send_rank);
  void SendReactors(std::vector<size_t> send_reactor_idxs, size_t recv_rank);
#endif

  void ProcessPerformance();

  void DumpReactor(std::string tag, int id, double T, double P,
                   double rc, double rg, double* mf);
};

#endif
