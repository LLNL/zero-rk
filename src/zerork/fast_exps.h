#ifndef ZERORK_FAST_EXPS_H
#define ZERORK_FAST_EXPS_H

namespace zerork {

extern const char* expType;
extern void (*fast_vec_exp)(double *, int);
extern double (*fast_exp)(double);

} // namespace zerork

#endif
