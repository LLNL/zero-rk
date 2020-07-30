#ifndef SET_INITIAL_CONDITIONS_H_
#define SET_INITIAL_CONDITIONS_H_

#include "flame_params.h"

// Set all grid points to the inlet conditions including temperature
void SetConstantInlet(FlameParams &flame_params, double *y);

#endif
