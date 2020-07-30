#ifndef SET_INITIAL_CONDITIONS_H_
#define SET_INITIAL_CONDITIONS_H_

#include "flame_params.h"

// Set initial conditions
void SetInitialCompositionAndWallTemp(FlameParams &flame_params, double *y, double *time);

#endif
