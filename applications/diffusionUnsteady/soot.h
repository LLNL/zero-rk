#ifndef SOOT_H_
#define SOOT_H_

#include "flame_params.h"
#include <nvector/nvector_parallel.h> // serial N_Vector types, fcts., and macros

void InitializePAH(void *user_data);

double CollisionEfficiency(const int i,
                           const int j,
                           const double temperature,
                           const double density,
                           const double molecular_mass[]);

void ComputeDimerProdRate(void *user_data,
                          const double state[],
                          double dimer_prod_rate[]);

double GetBetaDimer(const double temperature,
                    const double i,
                    const double j);

void CalcRhoDot(void *user_data,
                const double state[],
                double rhodot);

void UpdateProductionRates(void *user_data,
                           const double state[],
                           double prod_rate[]);

void ComputeDimerParticles(void *user_data,
		           const double state[]);

double GetBetaNucleation(const double temperature,
                         const double i);

void WriteDimerProdRate(void *user_data,
                        const double state[]);

void UpdateDimerProdRate(void *user_data,
                        const double state[]);

#endif
