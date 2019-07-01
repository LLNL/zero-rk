#ifndef ZERORK_FAST_EXPS_H
#define ZERORK_FAST_EXPS_H

/*
    Useful macro definitions for memory alignment:
        http://homepage1.nifty.com/herumi/prog/gcc-and-vc.html#MIE_ALIGN       
 */

#ifdef _MSC_VER
    #include <malloc.h>
#else
    #include <stdlib.h>
namespace zerork {
    static inline void *_aligned_malloc(size_t size, size_t alignment)
    {
        void *p;
        int ret = posix_memalign(&p, alignment, size);
        return (ret == 0) ? p : 0;
    }
} // namespace zerork
#endif

namespace zerork {

extern void (*fast_vec_exp)(double *, int);
extern double (*fast_exp)(double);

} // namespace zerork

#endif
