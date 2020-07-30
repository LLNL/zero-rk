#ifndef SEQUENTIAL_FILE_MATRIX_H_
#define SEQUENTIAL_FILE_MATRIX_H_

#include <string>

namespace zerork {
namespace utilities {

class SequentialFileMatrixRead
{
 public:
  SequentialFileMatrixRead(const char filename[]);
  ~SequentialFileMatrixRead();

  int GetNextRow(double data[]);

  // basic accessors
  int next_row_id() const {return next_row_id_;}
  int num_rows() const {return num_rows_;}
  int num_columns() const {return num_columns_;}
  int max_buffer_rows() const {return max_buffer_rows_;}
  const char * filename() const {return filename_.c_str();}

 private:
  int LoadBuffer();

  int next_row_id_;
  int num_rows_;
  int num_columns_;
  int max_buffer_rows_;
  double *next_buffer_ptr_;
  
  double *matrix_buffer_;

  std::string filename_;  

  FILE *fptr_;
};

} // end namespace utilities
} //end namespace zerork

#endif
