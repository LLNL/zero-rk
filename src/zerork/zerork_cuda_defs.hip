#include "zerork_cuda_defs.h"

#include <cstdio>
#include <cstdlib>

namespace zerork {

#ifdef EXIT_THROWS_EXCEPTION
  // create a local function to overide the system exit and throw an exception
  // with the status integer.
  static void exit(int status) {throw status;}
#endif // EXIT_THROWS_EXCEPTION

void checkCudaError(hipError_t err, const char * msg)
{
    if(!err)
        return;

    if(!msg)
    {
      std::printf("%s\n",hipGetErrorString(err));
    } else {
      std::printf("%s:  %s\n", msg, hipGetErrorString(err));
    }
#ifndef EXIT_THROWS_EXCEPTION
    std::exit(EXIT_FAILURE);
#else
    exit(EXIT_FAILURE);
#endif
    
}

}

