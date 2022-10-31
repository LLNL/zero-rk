#include <stdlib.h>
#include <stdio.h>
#include <algorithm> //std::min

#include "sequential_file_matrix.h"

namespace zerork {
namespace utilities {

SequentialFileMatrixRead::SequentialFileMatrixRead(const char filename[])
{
  int num_read;

  // open file
  fptr_ = fopen(filename,"r");
  if(fptr_ == NULL) {
    printf("ERROR: In SequentialFileMatrixRead(...),\n");
    printf("       could not open file named %s for read\n",
           filename);
    printf("       Exiting now.\n"); fflush(stdout);
    exit(-1);
  }
  filename_ = std::string(filename);

  // read the matrix size
  num_rows_    = -1;
  num_columns_ = -1;
  num_read = fscanf(fptr_,"%d%d",&num_rows_,&num_columns_);
  if(num_rows_ <= 0 || num_columns_ <= 0 || num_read != 2) {
    printf("ERROR: In SequentialFileMatrixRead(...),\n");
    printf("       failed to read the matrix size line.\n");
    printf("         Matrix file          = %s\n",filename_.c_str());
    printf("         No. of integers read = %d\n",num_read);
    printf("         No. of rows          = %d\n",num_rows_);
    printf("         No. of columns       = %d\n",num_columns_);
    printf("       Exiting now.\n"); fflush(stdout);
    exit(-1);
  }

  // compute number of buffer rows
  const size_t MAX_BUFFER_MEMORY = 1024*1024*10; //10 MB for buffer
  max_buffer_rows_ = MAX_BUFFER_MEMORY/(sizeof(double)*num_columns_);
  if(max_buffer_rows_ == 0) max_buffer_rows_ = 1;
  max_buffer_rows_ = std::min(max_buffer_rows_,num_rows_);

  // allocate buffer
  matrix_buffer_ = new double[max_buffer_rows_*num_columns_];
  if(matrix_buffer_ == NULL)  {
    printf("ERROR: In SequentialFileMatrixRead(...),\n");
    printf("       failed to allocate matrix buffer.\n");
    printf("         Matrix file = %s\n",filename_.c_str());
    printf("       Exiting now.\n"); fflush(stdout);
    exit(-1);    
  }

  // set the next global row index to be read
  next_row_id_ = 0;
  // load the first matrix buffer from the file
  if(int flag = LoadBuffer() != 0) {
    printf("ERROR: In SequentialFileMatrixRead(...),\n");
    printf("       failed to load matrix buffer.\n");
    printf("         Matrix file = %s\n",filename_.c_str());
    printf("         Error code  = %d\n",flag);
    printf("       Exiting now.\n"); fflush(stdout);
    exit(-1);
  }
}


SequentialFileMatrixRead::~SequentialFileMatrixRead()
{
  if(matrix_buffer_ != NULL) {
    delete [] matrix_buffer_;
  }

  if(fptr_ != NULL) {
    fclose(fptr_);
  }
}

// LoadBuffer()
//  pre-requisites:
//     private members set:  num_rows_
//                           num_columns_
//                           next_row_id_
//                           max_buffer_rows_
//                           fptr_
//                           filename_
//     private members allocated: matrix_buffer_
int SequentialFileMatrixRead::LoadBuffer()
{
  int num_read;
  const int row_limit = ((max_buffer_rows_ < num_rows_-next_row_id_) ?
                         max_buffer_rows_ : num_rows_-next_row_id_);
  const int column_limit = num_columns_;

  for(int j=0; j<row_limit; ++j) {
    for(int k=0; k<column_limit; ++k) {

      num_read = fscanf(fptr_,"%lf",&matrix_buffer_[j*column_limit+k]);

      if(num_read != 1) {
        printf("WARNING: In SequentialFileMatrixRead::LoadBuffer(),\n");
        printf("         failed to load buffer element.\n");
        printf("           Matrix file             = %s\n",filename_.c_str());
        printf("           Global row    (index-0) = %d\n",next_row_id_+j);
        printf("           Buffer row    (index-0) = %d\n",j);
        printf("           Buffer column (index-0) = %d\n",k);
        printf("         Exiting LoadBuffer().\n"); fflush(stdout);
        return -1;
      }      
    }
  }

  // reset the pointer to the next read location from the buffer
  next_buffer_ptr_ = &matrix_buffer_[0];
  return 0; // successfully loaded 
}

// GetNextRow(double data[])
//  pre-requisites:
//     private members set:  num_rows_
//                           num_columns_
//                           next_row_id_
//                           max_buffer_rows_
//                           fptr_
//                           filename_
//                           matrix_buffer_ <- LoadBuffer()
//                           next_buffer_ptr_
//
int SequentialFileMatrixRead::GetNextRow(double data[])
{
  const int column_limit = num_columns_;

  if(next_row_id_ >= num_rows_) {
    printf("WARNING: In SequentialFileMatrixRead::GetNextRow(...),\n");
    printf("         next row >= total rows.\n");
    printf("           Matrix file          = %s\n",filename_.c_str());
    printf("           next  row  = %d\n",next_row_id_);
    printf("           total rows = %d\n",num_rows_);
    printf("         Returning all zero data.\n");
    for(int j=0; j<column_limit; ++j) {
      data[j] = 0.0;
    }
    return -1;
  }

  for(int j=0; j<column_limit; ++j) {
    data[j] = next_buffer_ptr_[j];
  }
  ++next_row_id_;

  if(next_row_id_ < num_rows_) {
    // load a new buffer or advance the buffer pointer
    if(next_row_id_ % max_buffer_rows_ == 0) {
      if(int flag = LoadBuffer() != 0) {
        printf("ERROR: In SequentialFileMatrixRead::GetNextRow(...),\n");
        printf("       failed to load matrix buffer.\n");
        printf("         Matrix file = %s\n",filename_.c_str());
        printf("         Error code  = %d\n",flag);
        printf("       Exiting now.\n"); fflush(stdout);
        exit(-1);
      }
    } else {
      next_buffer_ptr_ += num_columns_;
    }
  }
  return next_row_id_;
}

} // end namespace utilities
} //end namespace zerork
