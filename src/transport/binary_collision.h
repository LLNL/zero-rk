#ifndef BINARY_COLLISION_H_
#define BINARY_COLLISION_H_

#include "utilities/sort_vector.h"

#include "zerork/mechanism.h"

namespace transport {

typedef struct {
  double epsilon_k;
  double sigma;      
} LJParameters;

void EstimateLJParameters(const double molecular_weight,
                          LJParameters *lj);
double GetHardSphereDiameter(const double viscosity_ref,
                             const double temperature_ref,
                             const double molecular_weight);
                             
enum LJMethod { INITIAL = -1, MW_ESTIMATE = 1, INTERNAL = 2, EXTERNAL = 3 };

class BinaryCollision
{
 public:
  BinaryCollision(zerork::mechanism *mech,
                  const double temp_ref);

  int num_binary_steps() const {return num_binary_steps_;}
  double temperature_ref() const {return temperature_ref_;} 
  int binary_step_index(const int id) const {return binary_step_index_[id];}

  double GetHSRateMultiplier(const int species_id1,
                             const int species_id2,
                             const double temperature);
  double GetStepProbability(const double temperature,
                            const double pressure,
                            double step_probability[]);
  void GetCollisionRateCoefficient(const double temperature,
                                   double step_collision_rate[]);

  void GetSpeciesReport(std::string *report);
  int GetStepProbabilityReport(const double pressure,
                               const double start_temperature,
                               const double end_temperature,
                               const int num_scan_points,
                               const double min_probability,
                               std::string *report);

 private:
  int BuildBinaryStepList();
  void SetEstimatedLJParameters();
  void SetInternalDatabase();
  void SetInternalLJParameters();
  void SetHSDiameter();
  void SetHSCollisionRateConstants();
  void WriteBinaryStepScanInBlocks(const int block_size,
                       const int num_scans,
                       const double scan_tempreature[],
		       const std::vector<zerork::utilities::SortableVector> *probability,
                       const double min_probability,
                       std::string *report);
  int num_binary_steps_;
  double temperature_ref_;

  zerork::mechanism *mechp_;

  std::vector<double> molecular_weight_;
  std::vector<double> viscosity_ref_;
  std::vector<double> hard_sphere_diameter_;
  std::vector<double> hard_sphere_collision_rate_constant_;
  std::vector<LJParameters> lennard_jones_;
  std::vector<LJMethod> method_;
  std::vector<int> binary_step_index_;

  std::map<std::string,LJParameters> internal_database_;

};
} // end of namespace transport
#endif
