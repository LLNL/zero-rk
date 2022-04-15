
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <iomanip>
#include <stdexcept>

#include "zerork_reactor_manager.h"
#include "ZeroRKCFDPluginIFP.h"

#include "solver_cvode.h"
#include "solver_seulex.h"
#include "reactor_constant_volume.h"
#include "reactor_constant_pressure.h"

#include "utility_funcs.h"

using zerork::getHighResolutionTime;

ZeroRKReactorManager::ZeroRKReactorManager(const char* input_filename,
                                           const char* mech_filename,
                                           const char* therm_filename)
  :
      ZeroRKReactorManagerBase(input_filename, mech_filename, therm_filename)
{
  n_calls_ = 0;
  n_cycle_ = 0;
  rank_ = 0;
  sorted_rank_ = 0;
  nranks_ = 1;
  root_rank_ = 0;
  load_balance_ = 1;
  dump_reactors_ = false;

#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&(rank_));
  MPI_Comm_size(MPI_COMM_WORLD,&(nranks_));
#endif

  //Default options
  //
  //Manager Options
  int_options_["verbosity"] = 4;
  int_options_["sort_reactors"]= 1;

  //Solver options
  int_options_["max_steps"] = 1000000;
  int_options_["dense"] = 0;
  int_options_["analytic"] = 1;
  int_options_["iterative"] = 1;
  int_options_["integrator"] = 0;
  int_options_["abstol_dens"] = 0;
  double_options_["rel_tol"] = 1.0e-8;
  double_options_["abs_tol"] = 1.0e-20;
  double_options_["eps_lin"] = 1.0e-3;
  double_options_["nonlinear_convergence_coeff"] = 0.05;
  double_options_["preconditioner_threshold"] = 1.0e-3;
  double_options_["max_dt"] = 0.05;


  //Reactor Options
  int_options_["constant_volume"] = 1;
  double_options_["reference_temperature"] = 1.0;
  double_options_["delta_temperature_ignition"] = 0.0;
  double_options_["min_mass_fraction"] = 1.0e-30;

  std::string reactor_timing_log_filename("/dev/null");

  //Overrides from input file
  ZeroRKCFDPluginIFP inputFileDB(input_filename);
  int_options_["verbosity"] = inputFileDB.verbosity();
  int_options_["sort_reactors"] = inputFileDB.sort_reactors();

  int_options_["max_steps"] = inputFileDB.max_steps();
  int_options_["dense"] = inputFileDB.dense();
  int_options_["analytic"] = inputFileDB.analytic();
  int_options_["iterative"] = inputFileDB.iterative();
  int_options_["integrator"] = inputFileDB.integrator();
  int_options_["abstol_dens"] = inputFileDB.abstol_dens();
  double_options_["abs_tol"] = inputFileDB.absolute_tolerance();
  double_options_["rel_tol"] = inputFileDB.relative_tolerance();
  double_options_["eps_lin"] = inputFileDB.eps_lin();
  double_options_["nonlinear_convergence_coeff"] = inputFileDB.nonlinear_convergence_coeff();
  double_options_["preconditioner_threshold"] = inputFileDB.preconditioner_threshold();
  double_options_["max_dt"] = inputFileDB.max_dt();

  int_options_["constant_volume"] = inputFileDB.constant_volume();
  double_options_["reference_temperature"] = inputFileDB.reference_temperature();
  double_options_["delta_temperature_ignition"] = inputFileDB.delta_temperature_ignition();
  double_options_["min_mass_fraction"] = inputFileDB.min_mass_fraction();

  n_reactors_max_ = inputFileDB.n_reactors_max();
  n_reactors_min_ = inputFileDB.n_reactors_min();
#ifdef USE_MPI
  load_balance_ = inputFileDB.load_balance();
  load_balance_noise_ = inputFileDB.load_balance_noise();
#endif
  if(nranks_ == 1) {
    load_balance_ = 0;
  }
  dump_reactors_ = (inputFileDB.dump_reactors() != 0);

  reactor_timing_log_filename = inputFileDB.reactor_timing_log();
  
  if(rank_ == 0) { 
    reactor_log_file_.open(reactor_timing_log_filename);
  }

  T_other_.clear();
  P_other_.clear();
  dpdt_other_.clear();
  e_src_other_.clear();
  y_src_other_.clear();
  mf_other_.clear();
  rc_other_.clear();
  rg_other_.clear();
  reactor_ids_other_.clear();
  root_times_other_.clear();

  dpdt_owned_ = true;
  rc_owned_ = true;
  rg_owned_ = true;
  root_times_owned_ = true;
  y_src_defined_ = false;
  e_src_defined_ = false;
  dpdt_self_ = nullptr;
  e_src_self_ = nullptr;
  y_src_self_ = nullptr;
  rc_self_ = nullptr;
  rg_self_ = nullptr;
  reactor_ids_defined_ = false;
  root_times_self_ = nullptr;

  cb_fn_ = nullptr;
  cb_fn_data_ = nullptr;

  const char* cklog_filename = inputFileDB.mechanism_parsing_log().c_str();
  if(rank_ != root_rank_) {
     cklog_filename = "/dev/null";
  }

  mech_ptr_ = std::make_shared<zerork::mechanism>(mech_filename, therm_filename, cklog_filename);
}

#ifndef USE_MPI
int ZeroRKReactorManager::RecvReactors(size_t send_rank) {
  return 0;
}
#else
int ZeroRKReactorManager::RecvReactors(size_t send_rank) {
  MPI_Status status;
  int n_recv_reactors;
  MPI_Recv(&n_recv_reactors,1, MPI_INT, send_rank, EXCHANGE_SEND_TAG_, MPI_COMM_WORLD,
           &status);
  if(n_recv_reactors == 0) return n_recv_reactors;

  size_t recv_idx = n_reactors_other_;
  size_t extra_tx_count = 0;
  n_reactors_other_ += n_recv_reactors;

  T_other_.resize(n_reactors_other_);
  P_other_.resize(n_reactors_other_);
  dpdt_other_.resize(n_reactors_other_);
  if(e_src_defined_) {
    e_src_other_.resize(n_reactors_other_);
    extra_tx_count += 1;
  }
  if(y_src_defined_) {
    y_src_other_.resize(n_reactors_other_*num_species_);
    extra_tx_count += num_species_;
  }
  rg_other_.resize(n_reactors_other_);
  rc_other_.resize(n_reactors_other_);
  root_times_other_.assign(n_reactors_other_, -1.0); //N.B. at start negative one by definition so we don't communicate it
  mf_other_.resize(n_reactors_other_*num_species_);

  std::vector<double> recv_buf(tx_count_per_reactor_+extra_tx_count);
  for(size_t i = 0; i < n_recv_reactors; ++i, ++recv_idx) {
     //Bring 'em in
     MPI_Recv(&recv_buf[0],tx_count_per_reactor_+extra_tx_count,
              MPI_DOUBLE, send_rank, EXCHANGE_SEND_TAG_,
              MPI_COMM_WORLD, &status);
     //Unpack 'em
     size_t buf_idx = 0;
     memcpy(&mf_other_[recv_idx*num_species_], &recv_buf[buf_idx],
            sizeof(double)*num_species_);
     buf_idx += num_species_;
     memcpy(&T_other_[recv_idx], &recv_buf[buf_idx], sizeof(double));
     buf_idx += 1;
     memcpy(&P_other_[recv_idx], &recv_buf[buf_idx], sizeof(double));
     buf_idx += 1;
     memcpy(&dpdt_other_[recv_idx], &recv_buf[buf_idx], sizeof(double));
     buf_idx += 1;
     if(e_src_defined_) {
       memcpy(&e_src_other_[recv_idx], &recv_buf[buf_idx], sizeof(double));
       buf_idx += 1;
     }
     if(y_src_defined_) {
       memcpy(&y_src_other_[recv_idx*num_species_], &recv_buf[buf_idx],
              sizeof(double)*num_species_);
       buf_idx += num_species_;
     }
     memcpy(&rg_other_[recv_idx], &recv_buf[buf_idx], sizeof(double));
     buf_idx += 1;
     memcpy(&rc_other_[recv_idx], &recv_buf[buf_idx], sizeof(double));
  }
  if(reactor_ids_defined_) {
     reactor_ids_other_.resize(n_reactors_other_);
     MPI_Recv(&reactor_ids_other_[0],n_reactors_other_,
              MPI_INT, send_rank, EXCHANGE_SEND_TAG_,
              MPI_COMM_WORLD, &status);
  }
  return n_recv_reactors;
}
#endif


#ifndef USE_MPI
void ZeroRKReactorManager::SendReactors(std::vector<size_t> send_reactor_idxs, size_t recv_rank) {}
#else
void ZeroRKReactorManager::SendReactors(std::vector<size_t> send_reactor_idxs, size_t recv_rank) {
  int send_nreactors = (int) send_reactor_idxs.size();
  MPI_Send(&send_nreactors,1,
           MPI_INT, (int) recv_rank, EXCHANGE_SEND_TAG_, MPI_COMM_WORLD);

  if(send_nreactors == 0) {
    return;
  }

  size_t extra_tx_count = 0;
  if(e_src_defined_) extra_tx_count+=1;
  if(y_src_defined_) extra_tx_count+=num_species_;
  std::vector<double> send_buf(tx_count_per_reactor_+extra_tx_count);
  for(int i = 0; i < send_nreactors; ++i) {
    int buf_idx = 0;
    int send_idx = send_reactor_idxs[i];
    memcpy(&send_buf[buf_idx], &mf_self_[send_idx*num_species_stride_],
           sizeof(double)*num_species_);
    buf_idx += num_species_;
    memcpy(&send_buf[buf_idx], &T_self_[send_idx], sizeof(double));
    buf_idx += 1;
    memcpy(&send_buf[buf_idx], &P_self_[send_idx], sizeof(double));
    buf_idx += 1;
    memcpy(&send_buf[buf_idx], &dpdt_self_[send_idx], sizeof(double));
    buf_idx += 1;
    if(e_src_defined_) {
      memcpy(&send_buf[buf_idx], &e_src_self_[send_idx], sizeof(double));
      buf_idx += 1;
    }
    if(y_src_defined_) {
      memcpy(&send_buf[buf_idx], &y_src_self_[send_idx*num_species_stride_],
             sizeof(double)*num_species_);
      buf_idx += num_species_;
    }
    memcpy(&send_buf[buf_idx], &rg_self_[send_idx], sizeof(double));
    buf_idx += 1;
    memcpy(&send_buf[buf_idx], &rc_self_[send_idx], sizeof(double));
    MPI_Send(&send_buf[0], tx_count_per_reactor_+extra_tx_count,
             MPI_DOUBLE, (int) recv_rank, EXCHANGE_SEND_TAG_, MPI_COMM_WORLD);
  }
  if(reactor_ids_defined_) {
     std::vector<int> send_ids(send_nreactors);
     for(int i = 0; i < send_nreactors; ++i) {
         int send_idx = send_reactor_idxs[i];
         send_ids[i] = reactor_ids_self_[send_idx];
     }
     MPI_Send(&send_ids[0], send_nreactors, MPI_INT,
              (int) recv_rank, EXCHANGE_SEND_TAG_, MPI_COMM_WORLD);
  }
}
#endif

void ZeroRKReactorManager::SetInputVariables(int n_cycle,
                                             double time,
                                             double dt,
                                             int n_reactors, 
                                             double* T,
                                             double* P,
                                             double* mf)

{
  n_cycle_ = n_cycle;
  dt_calc_ = dt;
  n_reactors_self_ = n_reactors;
  n_reactors_other_ = 0;
  T_self_ = T;
  P_self_ = P;
  mf_self_ = mf;

  if(dpdt_owned_) {
    dpdt_default_.resize(n_reactors, 0.0);
    dpdt_self_ = &dpdt_default_[0];
  }
  if(rg_owned_) {
    rg_default_.resize(n_reactors, 0.0);
    rg_self_ = &rg_default_[0];
  }
  if(rc_owned_) {
    rc_default_.resize(n_reactors, 1.0);
    rc_self_ = &rc_default_[0];
  }
  if(root_times_owned_) {
    root_times_default_.assign(n_reactors, -1.0);
    root_times_self_ = &root_times_default_[0];
  }

  sorted_reactor_idxs_.assign(n_reactors_self_,0);
  std::iota(sorted_reactor_idxs_.begin(), sorted_reactor_idxs_.end(), 0);
  if(int_options_["sort_reactors"]) {
      sorted_reactor_idxs_ = sort_indexes_pointer(n_reactors_self_, rc_self_);
  }
}

void ZeroRKReactorManager::SetAuxFieldPointer(zerork_field_type ft, double* field_pointer) {
  if(ft == ZERORK_FIELD_DPDT) {
      dpdt_self_ = field_pointer;
      dpdt_default_.clear();
      dpdt_owned_ = false;
  } else if (ft == ZERORK_FIELD_COST) {
      rc_self_ = field_pointer;
      rc_default_.clear();
      rc_owned_ = false;
  } else if (ft == ZERORK_FIELD_Y_SRC) {
      y_src_self_ = field_pointer;
      y_src_defined_ = true;
  } else if (ft == ZERORK_FIELD_E_SRC) {
      e_src_self_ = field_pointer;
      e_src_defined_ = true;
  } else if (ft == ZERORK_FIELD_IGNITION_TIME) {
      root_times_self_ = field_pointer;
      root_times_default_.clear();
      root_times_owned_ = false;
  } else {
      throw std::invalid_argument("Unknown zerork_field_type.");
  }
}

void ZeroRKReactorManager::SetReactorIDs(int* reactor_ids) {
  reactor_ids_self_ = reactor_ids;
  reactor_ids_defined_ = true;
}

void ZeroRKReactorManager::SetCallbackFunction(zerork_callback_fn fn, void* cb_fn_data) {
  cb_fn_ = fn;
  cb_fn_data_ = cb_fn_data;
}

void ZeroRKReactorManager::FinishInit() {
  if(n_calls_ == 0) {
    num_species_ = mech_ptr_->getNumSpecies();
    num_species_stride_ = num_species_;

    tx_count_per_reactor_ = 5 + num_species_;

#ifdef USE_MPI
    rank_weights_.assign(nranks_,0.0);
    double my_weight = 1.0;
    MPI_Allgather(&my_weight,1,MPI_DOUBLE,&rank_weights_[0],1,MPI_DOUBLE,MPI_COMM_WORLD);
    double sum_weights = 0;
    for(int i = 0 ; i < nranks_; ++i) {
      sum_weights += rank_weights_[i];
    }
    double weight_factor = nranks_/sum_weights;
    for(int i = 0 ; i < nranks_; ++i) {
      rank_weights_[i] *= weight_factor;
    }
#endif

    all_time_ranks_.assign(nranks_,0);
    n_reactors_solved_ranks_.assign(nranks_,0);
    n_weight_updates_ = 0;

    if(rank_ == 0) {
      if(int_options_["verbosity"] > 1) {
        printf("* %-25s%-31s *\n", "Zero-RK Lib Build Date: ",__DATE__);
      }

      // Timing log file
      reactor_log_file_ << "#";
      reactor_log_file_ << std::setw(11) << "solve_number";
      reactor_log_file_ << std::setw(17) << "reactors_solved";
      reactor_log_file_ << std::setw(17) << "n_steps_avg";
      reactor_log_file_ << std::setw(17) << "step_time";
      reactor_log_file_ << std::setw(17) << "avg_time";
      reactor_log_file_ << std::setw(17) << "max_time";
      reactor_log_file_ << std::endl;
      reactor_log_file_.flush();
    }

    num_species_stride_ = num_species_;
    if(int_options_.find("num_species_stride") != int_options_.end()) {
      num_species_stride_ = int_options_["num_species_stride"];
    }
  }
}

#ifndef USE_MPI
void ZeroRKReactorManager::LoadBalance()
{}; //pass
#else
void ZeroRKReactorManager::LoadBalance()
{
  if(nranks_ == 1 || !load_balance_) {
    return;
  }

  int n_weighted_reactors = 0;
  std::vector<int> weighted_reactors_on_rank(nranks_);
  for(int k = 0; k < n_reactors_self_; ++k) {
      if(rc_self_[k] <= 0.0) rc_self_[k] = 1.0;
      n_weighted_reactors += std::max((int)rc_self_[k],1)+load_balance_noise_;
  }

  MPI_Allgather(&n_weighted_reactors,1,MPI_INT,
                &weighted_reactors_on_rank[0],1,MPI_INT,MPI_COMM_WORLD);

  bool sending = false;
  bool recving = false;

  //Counts of how many reactors we want to add to ours to be balanced overall
  std::vector<int> reactors_wanted(nranks_);
  std::vector<int> reactor_deficit(nranks_);

  //Sparse matrix of reactors transmitted to/from ranks.
  //  (N.B. Column indicies are in reverse order from normal CSR)
  comm_mtx_row_sum_.assign(nranks_+1,0);
  comm_mtx_col_idx_.clear();
  comm_mtx_reactors_.clear();

  int total_weighted_reactors = 0;
  sorted_rank_ = -1;

  double sum_weights = 0.0;
  for(int i = 0; i < nranks_; ++i) {
    total_weighted_reactors += weighted_reactors_on_rank[i];
    sum_weights += rank_weights_[i];
  }

  int total_wanted = 0;
  double weight_factor = ((double)total_weighted_reactors)/sum_weights;
  for(int i = 0; i < nranks_; ++i) {
    reactors_wanted[i] = rank_weights_[i]*weight_factor;
    total_wanted += reactors_wanted[i];
  }

  for(int i = 0; i < nranks_; ++i) {
    reactor_deficit[i] = reactors_wanted[i] - weighted_reactors_on_rank[i];
    if(rank_==0 && int_options_["verbosity"] > 0) {
      printf("Rank[%d] reactors_wanted, weighted_reactors, reactor_deficit,"
            " rank_weight: %d, %d, %d, %f\n",
             i, reactors_wanted[i], weighted_reactors_on_rank[i],
             reactor_deficit[i],rank_weights_[i]);
    }
  }

  sorted_deficit_idxs_ = sort_indexes(reactor_deficit);
  for(int i = 0; i < nranks_; ++i) {
    //note where we are in the list
    if(sorted_deficit_idxs_[i] == rank_) {
      sorted_rank_ = i;
      break;
    }
  }
  assert(sorted_rank_ >= 0);

  //const int MIN_CHUNK = std::max(total_weighted_reactors/total_reactors,1);
  const int MIN_CHUNK = 1; //TODO
  int give_rank = nranks_-1; // smallest deficit (i.e. negative deficits)
                            // gets stolen from first
  for(int i = 0; i < nranks_; ++i) {
    comm_mtx_row_sum_[i+1] = comm_mtx_row_sum_[i]; //keep the running sum going
    double current_deficit = reactor_deficit[sorted_deficit_idxs_[i]];
    while(current_deficit >= MIN_CHUNK &&// while we stil want more
          give_rank > i) {                      //and there's more to be taken
      int extra_reactors_wanted    = current_deficit;
      int extra_reactors_available = -reactor_deficit[sorted_deficit_idxs_[give_rank]];
      if(extra_reactors_available >= MIN_CHUNK) {
        int reactors_taken = std::min(extra_reactors_available, extra_reactors_wanted);
        reactor_deficit[sorted_deficit_idxs_[i]] -= reactors_taken;
        reactor_deficit[sorted_deficit_idxs_[give_rank]] += reactors_taken;
        comm_mtx_row_sum_[i+1] += 1;
        comm_mtx_col_idx_.push_back(give_rank);
        comm_mtx_reactors_.push_back(reactors_taken);
      }
      if(-reactor_deficit[sorted_deficit_idxs_[give_rank]] < MIN_CHUNK) { //no more to give on this rank_
        give_rank -= 1;
      }
      if(reactor_deficit[sorted_deficit_idxs_[give_rank]] > 0) { //no more to give at all
        break;
      }
      current_deficit = reactor_deficit[sorted_deficit_idxs_[i]];
    }
  }
  for(int i = 0; i < comm_mtx_col_idx_.size(); ++i) {
    if(comm_mtx_col_idx_[i] == sorted_rank_) {
      sending = true;
      break;
    }
  }
  if((comm_mtx_row_sum_[sorted_rank_+1]-comm_mtx_row_sum_[sorted_rank_]) > 0) recving = true;

  //Exchange excess reactors
  if(recving) {
    for(int i = comm_mtx_row_sum_[sorted_rank_]; i < comm_mtx_row_sum_[sorted_rank_+1]; ++i) {
      int sorted_send_rank = comm_mtx_col_idx_[i];
      int send_rank = sorted_deficit_idxs_[sorted_send_rank]; //reactor_deficit[sorted_send_rank].first;
      int recv_nreactors = RecvReactors(send_rank);
      comm_mtx_reactors_[i] = recv_nreactors; //over-write weighted value with actual value
    }
  }

  if(sending) {
    int send_idx = n_reactors_self_;
    for(int ircv = 0; ircv < nranks_; ++ircv) {
      int recv_rank = sorted_deficit_idxs_[ircv]; //reactor_deficit[ircv].first;
      if(recv_rank == rank_) continue;
      for(int j = comm_mtx_row_sum_[ircv]; j < comm_mtx_row_sum_[ircv+1]; ++j) {
        int sorted_send_rank = comm_mtx_col_idx_[j];
        int send_rank = sorted_deficit_idxs_[sorted_send_rank]; //reactor_deficit[sorted_send_rank].first;
        if(send_rank == rank_) {
          int send_weighted_nreactors = comm_mtx_reactors_[j];
          int send_count = 0;
          std::vector<size_t> send_reactor_idxs;
          while(send_count < send_weighted_nreactors && send_idx>0) {
            send_idx -= 1;
            size_t sorted_reactor_idx = sorted_reactor_idxs_[send_idx];
            send_count += std::max((int)rc_self_[sorted_reactor_idx],1)+load_balance_noise_;
            send_reactor_idxs.push_back(sorted_reactor_idx);
          }
          size_t send_nreactors = send_reactor_idxs.size();
          comm_mtx_reactors_[j] = send_nreactors; //over-write weighted value
          n_reactors_other_ -= send_nreactors;
          SendReactors(send_reactor_idxs, recv_rank);
        }
      }
    }
  }

  if(n_reactors_other_ > 0) {
    if(int_options_["sort_reactors"]) {
      //To maintain sorting we need to sort again
      std::vector<double> rc_combined(rc_self_, rc_self_ + n_reactors_self_);
      rc_combined.insert(rc_combined.end(),rc_other_.begin(), rc_other_.end());
      sorted_reactor_idxs_ = sort_indexes(rc_combined);
    } else {
      //Need to resize to full length
      sorted_reactor_idxs_.assign(n_reactors_self_+n_reactors_other_,0);
      std::iota(sorted_reactor_idxs_.begin(), sorted_reactor_idxs_.end(), 0);
    }
  }
}
#endif //USE_MPI

void ZeroRKReactorManager::SolveReactors()
{


  n_calls_++;
  sum_reactor_time_ = 0.0;
  n_steps_ = 0;
  n_solve_ = 0;
  for(int j = 0; j < n_reactors_self_; ++j) {
     rg_self_[j] = 0.0;
  }
  if(n_reactors_other_ > 0) {
     rg_other_.assign(n_reactors_other_, 0.0);
  }

  int n_reactors_self_calc = n_reactors_self_ + n_reactors_other_;

  std::vector<double *> T_ptrs(n_reactors_self_calc);
  std::vector<double *> P_ptrs(n_reactors_self_calc);
  std::vector<double *> dpdt_ptrs(n_reactors_self_calc);
  std::vector<double *> e_src_ptrs(n_reactors_self_calc);
  std::vector<double *> y_src_ptrs(n_reactors_self_calc);
  std::vector<double *> rc_ptrs(n_reactors_self_calc);
  std::vector<double *> rg_ptrs(n_reactors_self_calc);
  std::vector<double *> mf_ptrs(n_reactors_self_calc);
  std::vector<double *> root_times_ptrs(n_reactors_self_calc);
  std::vector<int *> reactor_id_ptrs(n_reactors_self_calc);
  for(int j = 0; j < n_reactors_self_calc; ++j) {
    int j_sort = sorted_reactor_idxs_[j];
    if(j_sort < n_reactors_self_) {
      T_ptrs[j] = &T_self_[j_sort];
      P_ptrs[j] = &P_self_[j_sort];
      dpdt_ptrs[j] = &dpdt_self_[j_sort];
      if(e_src_defined_) {
        e_src_ptrs[j] = &e_src_self_[j_sort];
      }
      if(y_src_defined_) {
        y_src_ptrs[j] = &y_src_self_[j_sort*num_species_stride_];
      }
      if(reactor_ids_defined_) {
        reactor_id_ptrs[j] = &reactor_ids_self_[j_sort];
      }
      rc_ptrs[j] = &rc_self_[j_sort];
      rg_ptrs[j] = &rg_self_[j_sort];
      root_times_ptrs[j] = &root_times_self_[j_sort];
      mf_ptrs[j] = &mf_self_[j_sort*num_species_stride_];
    } else {
      j_sort -= n_reactors_self_;
      T_ptrs[j] = &T_other_[j_sort];
      P_ptrs[j] = &P_other_[j_sort];
      dpdt_ptrs[j] = &dpdt_other_[j_sort];
      if(e_src_defined_) {
        e_src_ptrs[j] = &e_src_other_[j_sort];
      }
      if(y_src_defined_) {
        y_src_ptrs[j] = &y_src_other_[j_sort*num_species_];
      }
      if(reactor_ids_defined_) {
        reactor_id_ptrs[j] = &reactor_ids_other_[j_sort];
      }
      rc_ptrs[j] = &rc_other_[j_sort];
      rg_ptrs[j] = &rg_other_[j_sort];
      root_times_ptrs[j] = &root_times_other_[j_sort];
      mf_ptrs[j] = &mf_other_[j_sort*num_species_];
    }
    if(dump_reactors_) {
      DumpReactor("pre", j, *T_ptrs[j], *P_ptrs[j], *rc_ptrs[j], mf_ptrs[j]);
    }
  }

  //Instantiate reactors on first call, after options are set
  if(!reactor_ptr_) {
    if(int_options_["constant_volume"] == 1) {
      reactor_ptr_ = std::make_unique<ReactorConstantVolume>(mech_ptr_);
    } else {
      reactor_ptr_ = std::make_unique<ReactorConstantPressure>(mech_ptr_);
    }
  }

  reactor_ptr_->SetIntOptions(int_options_);
  reactor_ptr_->SetDoubleOptions(double_options_);

  std::unique_ptr<SolverBase> solver;
  if(int_options_["integrator"] == 0) {
    solver.reset(new CvodeSolver(*reactor_ptr_));
  } else {
    solver.reset(new SeulexSolver(*reactor_ptr_));
  }
  solver->SetIntOptions(int_options_);
  solver->SetDoubleOptions(double_options_);
  if(cb_fn_ != nullptr && load_balance_ == 0 && n_reactors_self_calc == 1) {
    solver->SetCallbackFunction(cb_fn_, cb_fn_data_);
  }

  reactor_ptr_->SetIntOption("iterative",solver->Iterative());

  for(int k = 0; k < n_reactors_self_calc; ++k)
  {
    if(*rg_ptrs[k] == 0) {
      double start_time = getHighResolutionTime();
      double e_src_reactor = 0.0;
      if(e_src_defined_) {
        e_src_reactor = *e_src_ptrs[k];
      }
      double* y_src_reactor = nullptr;
      if(y_src_defined_) {
        y_src_reactor = y_src_ptrs[k];
      }
      int reactor_id = k;
      if(reactor_ids_defined_) {
        reactor_id = *reactor_id_ptrs[k];
      }
      reactor_ptr_->SetID(reactor_id);
      reactor_ptr_->InitializeState(0.0, 1, T_ptrs[k], P_ptrs[k],
                                    mf_ptrs[k], dpdt_ptrs[k],
                                    &e_src_reactor,
                                    y_src_reactor);
      int nsteps = solver->Integrate(dt_calc_);
      reactor_ptr_->GetState(T_ptrs[k], P_ptrs[k], mf_ptrs[k]);
      *root_times_ptrs[k] = reactor_ptr_->GetRootTime();
      double reactor_time = getHighResolutionTime() - start_time;
      *rc_ptrs[k] = nsteps;
      n_steps_ += nsteps;
      sum_reactor_time_ += reactor_time;
      ++n_solve_;
      if(dump_reactors_) {
        DumpReactor("postc", k, *T_ptrs[k], *P_ptrs[k],
                    *rc_ptrs[k], mf_ptrs[k]);
      }
    }
  }
}

#ifndef USE_MPI
void ZeroRKReactorManager::RedistributeResults()
{}
#else
void ZeroRKReactorManager::RedistributeResults()
{
  if(nranks_ == 1 || !load_balance_) {
    return;
  }
  //startTime = getHighResolutionTime();
  MPI_Barrier(MPI_COMM_WORLD);
  //double synchTime = getHighResolutionTime() - startTime;

  if(n_reactors_other_ > 0) {
    int send_idx = 0;
    for(int i = comm_mtx_row_sum_[sorted_rank_]; i < comm_mtx_row_sum_[sorted_rank_+1]; ++i)
    {
      int sorted_recv_rank = comm_mtx_col_idx_[i];
      int recv_rank = sorted_deficit_idxs_[sorted_recv_rank];
      int send_nreactors = comm_mtx_reactors_[i];
      if(send_nreactors > 0)
      {
        //TODO: send single buffer
        for(int k = 0; k < send_nreactors; ++k) {
          MPI_Send(&rg_other_[send_idx], 1, MPI_DOUBLE,
                   recv_rank, EXCHANGE_RETURN_TAG_,MPI_COMM_WORLD);
          MPI_Send(&rc_other_[send_idx], 1, MPI_DOUBLE,
                   recv_rank, EXCHANGE_RETURN_TAG_,MPI_COMM_WORLD);
          MPI_Send(&T_other_[send_idx], 1, MPI_DOUBLE,
                   recv_rank, EXCHANGE_RETURN_TAG_,MPI_COMM_WORLD);
          MPI_Send(&P_other_[send_idx], 1, MPI_DOUBLE,
                   recv_rank, EXCHANGE_RETURN_TAG_,MPI_COMM_WORLD);
          MPI_Send(&root_times_other_[send_idx], 1, MPI_DOUBLE,
                   recv_rank, EXCHANGE_RETURN_TAG_,MPI_COMM_WORLD);
          MPI_Send(&mf_other_[send_idx*num_species_], num_species_,
                   MPI_DOUBLE, recv_rank, EXCHANGE_RETURN_TAG_,MPI_COMM_WORLD);
          send_idx += 1;
        }
      }
    }
  }

  if(n_reactors_other_ < 0) {
    MPI_Status status;
    int recv_idx = n_reactors_self_;
    for(int ircv = 0; ircv < nranks_; ++ircv) {
      int send_rank = sorted_deficit_idxs_[ircv];
      if(send_rank == rank_) continue;
      for(int j = comm_mtx_row_sum_[ircv]; j < comm_mtx_row_sum_[ircv+1]; ++j) {
        int sorted_recv_rank = comm_mtx_col_idx_[j];
        int recv_rank = sorted_deficit_idxs_[sorted_recv_rank];
        if(recv_rank == rank_) {
          int recv_nreactors = comm_mtx_reactors_[j];
          for(int k = 0; k < recv_nreactors; ++k) {
            recv_idx -= 1;
            size_t sorted_recv_idx = sorted_reactor_idxs_[recv_idx];
            MPI_Recv(&rg_self_[sorted_recv_idx], 1, MPI_DOUBLE,
                     send_rank, EXCHANGE_RETURN_TAG_, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&rc_self_[sorted_recv_idx], 1, MPI_DOUBLE,
                     send_rank, EXCHANGE_RETURN_TAG_, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&T_self_[sorted_recv_idx], 1, MPI_DOUBLE,
                     send_rank, EXCHANGE_RETURN_TAG_, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&P_self_[sorted_recv_idx], 1, MPI_DOUBLE,
                     send_rank, EXCHANGE_RETURN_TAG_, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&root_times_self_[sorted_recv_idx], 1, MPI_DOUBLE,
                     send_rank, EXCHANGE_RETURN_TAG_, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&mf_self_[sorted_recv_idx*num_species_stride_], num_species_,
                     MPI_DOUBLE, send_rank, EXCHANGE_RETURN_TAG_, MPI_COMM_WORLD,
                     &status);
          }
        }
      }
    }
  }
  //commTime += getHighResolutionTime() - startTime;
}
#endif

void ZeroRKReactorManager::PostSolve() {
  UpdateRankWeights();
  ProcessPerformance();
}

void ZeroRKReactorManager::ProcessPerformance()
{
  //N.B. These are set in case of no-mpi or mpi with single rank
  n_reactors_solved_ranks_[rank_] = n_solve_;
  all_time_ranks_[rank_] = sum_reactor_time_;
#ifdef USE_MPI
  if(nranks_ > 1) {
    MPI_Gather(&n_solve_,1,MPI_INT,&n_reactors_solved_ranks_[0],1,MPI_INT,root_rank_,MPI_COMM_WORLD);
    MPI_Gather(&sum_reactor_time_,1,MPI_DOUBLE,&all_time_ranks_[0],1,MPI_DOUBLE,root_rank_,MPI_COMM_WORLD);

    int ri; //reduced int
    MPI_Reduce(&n_steps_,&ri,1,MPI_INT,MPI_SUM,root_rank_,MPI_COMM_WORLD);
    if(rank_ == 0) n_steps_ = ri;
  }
#endif
  if(rank_ == root_rank_) {
    double max_time = 0.0;
    double total_time = 0.0;
    int total_solved = 0;
    for(int i = 0; i < nranks_; ++i) {
       max_time = std::max(max_time,all_time_ranks_[i]);
       total_time += all_time_ranks_[i];
       total_solved += n_reactors_solved_ranks_[i];
    }
    double avg_time = total_time/nranks_;
    double nstep_avg = total_solved > 0 ? n_steps_/total_solved : 0;
    double per_step_time = n_steps_ > 0 ? total_time/n_steps_ : 0;

    //Print stats to file
    //reactor_log_file_ << std::setprecision(16);
    reactor_log_file_ << std::setw(13) <<  n_cycle_;
    reactor_log_file_ << std::setw(17) <<  total_solved;
    reactor_log_file_ << std::setw(17) <<  nstep_avg;
    reactor_log_file_ << std::setw(17) <<  per_step_time;
    reactor_log_file_ << std::setw(17) <<  avg_time;
    reactor_log_file_ << std::setw(17) <<  max_time;
    reactor_log_file_ << std::endl;
    reactor_log_file_.flush();

    if(int_options_["verbosity"] > 0) {
      double wasted_time = max_time - avg_time;
      for(int i = 0; i < nranks_; ++i) {
        printf("Rank %d calculated %d reactors in %f seconds.\n",i,n_reactors_solved_ranks_[i],all_time_ranks_[i]);
      }
      printf("Max Time, Avg Time, Wasted Time = %f, %f, %f\n",max_time,avg_time,wasted_time);
    }
  }
  //commTime += getHighResolutionTime() - startTime;
}



#ifndef USE_MPI
void ZeroRKReactorManager::UpdateRankWeights() {}
#else
void ZeroRKReactorManager::UpdateRankWeights() {
  if(nranks_ == 1 || !load_balance_) {
    return;
  }
  //startTime = getHighResolutionTime();
  //Update rank_weights_
  if(rank_==0) {
    double total_all_time_ranks = 0.0;
    for(int i = 0; i < nranks_; ++i) {
      total_all_time_ranks += all_time_ranks_[i];
    }
    double avg_time_ranks = total_all_time_ranks/nranks_;
    const double alpha = 0.5; //Don't over correct
    double sum_weights = 0.0;
    for(int i = 0; i < nranks_; ++i) {
      double last_weight = rank_weights_[i];
      double new_weight;
      if(n_reactors_solved_ranks_[i] > 0) {
          new_weight = avg_time_ranks/all_time_ranks_[i]*last_weight;
      } else {
        new_weight = last_weight;
      }
      rank_weights_[i] =(1.0-alpha)*last_weight + alpha*new_weight;
      sum_weights += rank_weights_[i];
    }
    double weight_factor = nranks_/sum_weights;
    for(int i = 0 ; i < nranks_; ++i) {
      rank_weights_[i] *= weight_factor;
      if(int_options_["verbosity"] >= 2) {
        printf("RANK[%d] weight: %f\n",i,rank_weights_[i]);
      }
    }
  }

  MPI_Bcast(&rank_weights_[0],nranks_,MPI_DOUBLE,root_rank_,MPI_COMM_WORLD);
  n_weight_updates_ += 1;
  //commTime += getHighResolutionTime() - startTime;
}
#endif

void ZeroRKReactorManager::DumpReactor(std::string tag, int id, double T, double P,
                                       double rc, double* mf) {
      std::ofstream dump_file;
      std::ostringstream dump_file_name;
      dump_file_name << std::setfill('0') << "dumpfile_" << tag
                     << "_" << std::setw(6)
                     << n_calls_ << "_" << std::setw(6)
                     << rank_;// << "_" << std::setw(6) << id;
      dump_file.open(dump_file_name.str(), std::ofstream::out | std::ofstream::app);
      dump_file << id << std::endl;
      dump_file << std::setprecision(20);
      for(int k = 0 ; k < num_species_; ++k) {
        dump_file << mf[k] << std::endl;
      }
      dump_file << T << std::endl;
      dump_file << P << std::endl;
      dump_file << rc << std::endl;
      dump_file.close();
}
