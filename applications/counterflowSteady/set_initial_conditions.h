#ifndef SET_INITIAL_CONDITIONS_H_
#define SET_INITIAL_CONDITIONS_H_

#include "flame_params.h"
#include <nvector/nvector_parallel.h> // serial N_Vector types, fcts., and macros

// Set initial composition
//void SetInitialComposition(FlameParams &flame_params, double *y, double *ydot, double *id, double *time);

void SetInitialComposition(FlameParams &flame_params, N_Vector yvec, double *time);

#endif
