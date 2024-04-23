#ifndef GPU_ERR_CHECK_H
#define GPU_ERR_CHECK_H

#include <stdio.h>

//https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
//(accessed 20171220)
#define gpu_err_check(ans) { gpu_assert((ans), __FILE__, __LINE__); }
inline void gpu_assert(hipError_t code, const char *file, int line, bool abort=true)
{
   if (code != hipSuccess) 
   {
      fprintf(stderr,"gpu_assert: %s %s %d\n", hipGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#endif
