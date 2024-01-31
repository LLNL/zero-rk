
#include "mpi_utilities.h"

#include <stdio.h>

namespace zerork 
{
namespace utilities
{

mpi_counter::mpi_counter(MPI_Comm comm, int hostrank) {
    comm_ = comm;
    host_rank_ = hostrank;
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    if (rank_ == host_rank_) {
        data_.assign(size_,0);
        MPI_Win_create(&data_[0], size_*sizeof(int), sizeof(int),
                       MPI_INFO_NULL, comm_, &win_);
    } else {
        data_.clear();
        MPI_Win_create(nullptr, 0, 1,
                       MPI_INFO_NULL, comm_, &win_);
    }
    my_val_ = 0;
}

mpi_counter::~mpi_counter() {
  MPI_Win_free(&win_);
}

int mpi_counter::increment(int amount) {
    std::vector<int> vals(size_, 0);

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, host_rank_, 0, win_);

    for (int i=0; i<size_; i++) {
      if (i == rank_) {
        MPI_Accumulate(&amount, 1, MPI_INT, host_rank_, i, 1, MPI_INT, MPI_SUM,
                       win_);
      } else {
        MPI_Get(&vals[i], 1, MPI_INT, host_rank_, i, 1, MPI_INT, win_);
      }
    }

    MPI_Win_unlock(host_rank_, win_);

    vals[rank_] = my_val_;

    int val = 0;
    for (int i=0; i<size_; i++) {
      val += vals[i];
    }

    my_val_ += amount;

    return val; //value prior to increment
}


void mpi_counter::print() {
  if (rank_ == host_rank_) {
    for (int i=0; i<size_; i++) {
        printf("%2d ", data_[i]);
    }
    puts("");
  }
  fflush(stdout);
}

}
}
