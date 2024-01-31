#ifdef _WIN32
#define _USE_MATH_DEFINES //for M_PI
#endif
#include "math.h"

#include <zerork/constants.h>
#include "collision_integrals.h"

#include "binary_collision.h"

namespace transport {

// Use the Wang & Frenklach molecular weight correlations for PAH's as
// an estimate of the Lennard Jones 6-12 potential parameters. 
//
// Reference: H. Wang & M. Frenklach, "Transport Properties of Polycyclic
//            Aromatic Hydrocarbons for Flame Modeling," Combust. Flame,
//            96:163-170, 1994.
void EstimateLJParameters(const double molecular_weight,
                          LJParameters *lj)
{
  // epsilon_k is the well-depth divided by the Boltzmann constant and
  // has units of Kelvin, and the molecular weight is given in [g/mol] or
  // [kg/kmol]
  lj->epsilon_k = 37.15*pow(molecular_weight,0.58);
  // sigma is the lennard-jones collision diameter and corresponds to the 
  // distance when the inter-particle potential is zero.  Here sigma has units
  // of meters and  the molecular weight is given in [g/mol] or [kg/kmol]
  lj->sigma = 1.236e-10*pow(molecular_weight,1.0/3.0);
}

// Return the hard sphere diameter (in [m]) corresponding to the given
// inputs:
//   viscosity in [Pa-s]
//   temperature in [K] for the viscosity
//   molecular weight in [g/mol] or [kg/kmol]
//
// Eq 3.59 used from the book by Bird (Molecular Gas Dynamics and the Direct
// Simulation of Gas Flows, Clarendon Press, Oxford, 1994.) 
double GetHardSphereDiameter(const double viscosity_ref,
                             const double temperature_ref,
                             const double molecular_weight)
{
  return sqrt(5.0/(16.0*viscosity_ref)*sqrt(molecular_weight*
                                            zerork::KBoltzmann*
                                            temperature_ref/
                                            (M_PI*zerork::KAvogadroNumber)));
}

BinaryCollision::BinaryCollision(zerork::mechanism *mech,
                                 const double temp_ref)
{
  LJParameters initial_lj;

  mechp_ = mech;
  temperature_ref_ = temp_ref;

  molecular_weight_.assign(mechp_->getNumSpecies(),0.0);
  mechp_->getMolWtSpc(molecular_weight_.data());
  hard_sphere_diameter_.assign(mechp_->getNumSpecies(),0.0);
  hard_sphere_collision_rate_constant_.assign(mechp_->getNumSteps(),0.0);
  viscosity_ref_.assign(mechp_->getNumSpecies(),0.0);
  method_.assign(mechp_->getNumSpecies(),INITIAL);

  // initialize the Lennard-Jones parameter array
  // this should be over-written
  initial_lj.epsilon_k = 71.4;      // [K] values from R. Svehla NASA TR R-132
  initial_lj.sigma     = 3.798e-10; // [m] from 1962 for diatomic nitrogen
  lennard_jones_.assign(mechp_->getNumSpecies(),initial_lj);

  num_binary_steps_ = BuildBinaryStepList();

  SetEstimatedLJParameters();
  SetInternalLJParameters();

  // computes the reference viscosity and associated hard-sphere diameter
  SetHSDiameter();
  SetHSCollisionRateConstants();
}
// Returns the hard sphere collision rate multiplier Kmult, defined such that
// No. of collisions per second per m^3 = Kmult*n[1]*n[2]
// where n[1] and n[2] are the number density [particles/m^3] of species 1
// and 2 given by the species_id1 and species_id2.
//
// The hard sphere collision rate is given in eq 3.134 in T.I. Gombosi,
// Gaskinetic Theory, Cambridge University Press, 1994.               
double BinaryCollision::GetHSRateMultiplier(const int species_id1,
                                            const int species_id2,
                                            const double temperature)
{
  double Kmult = 1.0;
  double diameter_avg = 0.5*(hard_sphere_diameter_[species_id1]+
                             hard_sphere_diameter_[species_id2]);

  double reduced_mass = (molecular_weight_[species_id1]*
                         molecular_weight_[species_id2])/
                        (molecular_weight_[species_id1]+
		         molecular_weight_[species_id2])/
                        zerork::KAvogadroNumber;         
  // reduce Kmult by 1/2 to avoid double counting collisions when the particles
  // are the same
  if(species_id1 == species_id2) {
    Kmult = 0.5;
  }
  Kmult*=(M_PI*diameter_avg*diameter_avg); // multiply by the cross-section
  Kmult*=sqrt(8.0*zerork::KBoltzmann*temperature/(M_PI*reduced_mass));
  return Kmult;
}
double BinaryCollision::GetStepProbability(const double temperature,
                                           const double pressure,
                                           double step_probability[])
{
  int species_id1,species_id2, reaction_id;
  double step_rate_constant;
  double max_probability = 0.0;
  double total_concentration = pressure/(mechp_->getGasConstant()*temperature);
  std::vector<double> concentration;
  std::vector<double> Kfwd;
  std::vector<double> Krev;

  // The concentration vector is set to the average species concentration.
  // The actual species concentration should not affect the binary step rates.
  // The total concentration via the pressure will impact binary step rates
  // for PLOG reaction types.
  concentration.assign(mechp_->getNumSpecies(),
                       total_concentration/
                       static_cast<double>(mechp_->getNumSpecies()));
  Kfwd.assign(mechp_->getNumReactions(),0.0);
  Krev.assign(mechp_->getNumReactions(),0.0);
  
  mechp_->getKrxnFromTC(temperature,
                        concentration.data(),
                        Kfwd.data(),
                        Krev.data());

  for(int j=0; j<num_binary_steps_; ++j) { 

    species_id1 = mechp_->getSpecIdxOfStepReactant(binary_step_index_[j],0);
    species_id2 = mechp_->getSpecIdxOfStepReactant(binary_step_index_[j],1);

    reaction_id = mechp_->getRxnIdxOfStep(binary_step_index_[j]);

    if(binary_step_index_[j] == mechp_->getStepIdxOfRxn(reaction_id,1)) {
      step_rate_constant = Kfwd[reaction_id];
    } else {
      step_rate_constant = Krev[reaction_id];
    }
    // binary rate constant units [m^3/kmol/s]
    // convert to [m^3/particles/s]
    step_rate_constant *= (1.0/zerork::KAvogadroNumber);

    step_probability[j] = step_rate_constant/GetHSRateMultiplier(species_id1,
                                                                 species_id2,
                                                                 temperature);
    if(step_probability[j] > max_probability) {
      max_probability = step_probability[j];
    }
  }
  return max_probability;
}

void BinaryCollision::GetCollisionRateCoefficient(const double temperature,
                                                 double step_collision_rate[])
{
  const int num_steps = mechp_->getNumSteps();
  const double sqrt_temperature = sqrt(temperature);

  for(int j=0; j<num_steps; ++j) {
    step_collision_rate[j] = 
      sqrt_temperature*hard_sphere_collision_rate_constant_[j];
  }
}

int BinaryCollision::BuildBinaryStepList()
{
  int reaction_index;

  binary_step_index_.clear();
  for(int j=0; j<mechp_->getNumSteps(); ++j) {
    reaction_index = mechp_->getRxnIdxOfStep(j);

    if(mechp_->isFalloffReaction(reaction_index) == 0 &&
       mechp_->isThirdBodyReaction(reaction_index) == 0) {
 
      if(mechp_->getOrderOfStep(j) == 2) {
        binary_step_index_.push_back(j);
      }
    }
  }
  //printf("# Found %d binary steps\n",binary_step_list_.size());
  return binary_step_index_.size();
}

void BinaryCollision::SetEstimatedLJParameters()
{
  for(int j=0; j<mechp_->getNumSpecies(); ++j) {
    EstimateLJParameters(molecular_weight_[j],
                         &lennard_jones_[j]);
    method_[j] = MW_ESTIMATE;
  }
}
void BinaryCollision::SetInternalLJParameters()
{
  std::string lower_case_name;
  std::map<std::string,LJParameters>::iterator db_iter;

  SetInternalDatabase();
  for(int j=0; j<mechp_->getNumSpecies(); ++j) {
    lower_case_name = std::string(mechp_->getSpeciesName(j));
    for(unsigned int k=0; k<lower_case_name.length(); ++k) {
      char c = lower_case_name[k];
      if('A' <= c && c <= 'Z') {
        c+=('a'-'A');
        lower_case_name[k] = c;
      }
    } // now every upper case letter should be changed to lower case
    db_iter=internal_database_.find(lower_case_name);
    if(db_iter != internal_database_.end()) {
      lennard_jones_[j].epsilon_k = db_iter->second.epsilon_k;
      lennard_jones_[j].sigma = db_iter->second.sigma;
      method_[j] = INTERNAL;
    }
  }

}

void BinaryCollision::SetInternalDatabase()
{
  std::pair<std::string,LJParameters> add_lj;

  // Species from R. Svehla, NASA Tech Report TR R-132, 1962
  add_lj.first = "he";
  add_lj.second.sigma = 2.551e-10;
  add_lj.second.epsilon_k = 10.22;
  internal_database_.insert(add_lj);
  add_lj.first = "ar";
  add_lj.second.sigma = 3.542e-10;
  add_lj.second.epsilon_k = 93.3;
  internal_database_.insert(add_lj);
  add_lj.first = "h";
  add_lj.second.sigma = 2.708e-10;
  add_lj.second.epsilon_k = 37.0;
  internal_database_.insert(add_lj);
  add_lj.first = "h2";
  add_lj.second.sigma = 2.827e-10;
  add_lj.second.epsilon_k = 59.7;
  internal_database_.insert(add_lj);
  add_lj.first = "o";
  add_lj.second.sigma = 3.050e-10;
  add_lj.second.epsilon_k = 106.7;
  internal_database_.insert(add_lj);
  add_lj.first = "o2";
  add_lj.second.sigma = 2.467e-10;
  add_lj.second.epsilon_k = 106.7;
  internal_database_.insert(add_lj);
  add_lj.first = "oh";
  add_lj.second.sigma = 3.147e-10;
  add_lj.second.epsilon_k = 79.8;
  internal_database_.insert(add_lj);
  add_lj.first = "h2o";
  add_lj.second.sigma = 2.641e-10;
  add_lj.second.epsilon_k = 809.1;
  internal_database_.insert(add_lj);
  add_lj.first = "h2o2";
  add_lj.second.sigma = 4.196e-10;
  add_lj.second.epsilon_k = 289.3;
  internal_database_.insert(add_lj);
  add_lj.first = "n2";
  add_lj.second.sigma = 3.798e-10;
  add_lj.second.epsilon_k = 71.4;
  internal_database_.insert(add_lj);
  add_lj.first = "no";
  add_lj.second.sigma = 3.492e-10;
  add_lj.second.epsilon_k = 116.7;
  internal_database_.insert(add_lj);
  add_lj.first = "n2o";
  add_lj.second.sigma = 3.828e-10;
  add_lj.second.epsilon_k = 232.4;
  internal_database_.insert(add_lj);
  add_lj.first = "nh";
  add_lj.second.sigma = 3.312e-10;
  add_lj.second.epsilon_k = 65.3;
  internal_database_.insert(add_lj);
  add_lj.first = "nh3";
  add_lj.second.sigma = 2.900e-10;
  add_lj.second.epsilon_k = 558.3;
  internal_database_.insert(add_lj); 
  add_lj.first = "co"; // based on experiment
  add_lj.second.sigma = 3.690e-10;
  add_lj.second.epsilon_k = 91.7;
  internal_database_.insert(add_lj); 
  add_lj.first = "co2"; // based on experiment
  add_lj.second.sigma = 3.941e-10;
  add_lj.second.epsilon_k = 195.2;
  internal_database_.insert(add_lj); 
  add_lj.first = "co"; // based on experiment
  add_lj.second.sigma = 3.690e-10;
  add_lj.second.epsilon_k = 91.7;
  internal_database_.insert(add_lj); 
  add_lj.first = "ch4"; // based on experiment
  add_lj.second.sigma = 3.758e-10;
  add_lj.second.epsilon_k = 148.6;
  internal_database_.insert(add_lj); 
  add_lj.first = "c2h2"; // based on experiment
  add_lj.second.sigma = 4.033e-10;
  add_lj.second.epsilon_k = 231.8;
  internal_database_.insert(add_lj); 
  add_lj.first = "c2h4"; // based on experiment
  add_lj.second.sigma = 4.163e-10;
  add_lj.second.epsilon_k = 224.7;
  internal_database_.insert(add_lj); 
  add_lj.first = "c2h6"; // based on experiment
  add_lj.second.sigma = 4.443e-10;
  add_lj.second.epsilon_k = 215.7;
  internal_database_.insert(add_lj); 
  add_lj.first = "c3h8"; // based on experiment
  add_lj.second.sigma = 5.118e-10;
  add_lj.second.epsilon_k = 237.1;
  internal_database_.insert(add_lj); 
  add_lj.first = "nc4h10"; // based on experiment
  add_lj.second.sigma = 4.687e-10;
  add_lj.second.epsilon_k = 531.4;
  internal_database_.insert(add_lj); 
  add_lj.first = "ic4h10"; // based on experiment
  add_lj.second.sigma = 5.278e-10;
  add_lj.second.epsilon_k = 330.1;
  internal_database_.insert(add_lj); 
  add_lj.first = "nc5h12"; // based on experiment
  add_lj.second.sigma = 5.784e-10;
  add_lj.second.epsilon_k = 341.1;
  internal_database_.insert(add_lj); 
  add_lj.first = "nc6h14"; // based on experiment
  add_lj.second.sigma = 5.949e-10;
  add_lj.second.epsilon_k = 399.3;
  internal_database_.insert(add_lj); 
  add_lj.first = "c6h6"; // based on experiment
  add_lj.second.sigma = 5.349e-10;
  add_lj.second.epsilon_k = 412.3;
  internal_database_.insert(add_lj); 
 }
void BinaryCollision::SetHSDiameter()
{
  for(int j=0; j<mechp_->getNumSpecies(); ++j) {
    viscosity_ref_[j] = getViscosity_LJ(temperature_ref_,
                                        molecular_weight_[j],
                                        lennard_jones_[j].epsilon_k,
                                        lennard_jones_[j].sigma);
    hard_sphere_diameter_[j] = GetHardSphereDiameter(viscosity_ref_[j],
                                                 temperature_ref_,
                                                 molecular_weight_[j]);
  }
}
void BinaryCollision::SetHSCollisionRateConstants()
{
  const int num_binary_steps = (int)binary_step_index_.size();
  for(int j=0; j<num_binary_steps; ++j) {

    int step_id = binary_step_index_[j];
    int species0_id = mechp_->getSpecIdxOfStepReactant(step_id,0);
    int species1_id = mechp_->getSpecIdxOfStepReactant(step_id,1);
    
    double diameter = 0.5*(hard_sphere_diameter_[species0_id]+
                           hard_sphere_diameter_[species1_id]); // [m]
    double mass = molecular_weight_[species0_id]*
      molecular_weight_[species1_id]/
      (molecular_weight_[species0_id]+molecular_weight_[species1_id])/
      zerork::KAvogadroNumber; // [kg] reduced mass

    double collision_constant = 1.0;
    if(species0_id == species1_id) {
      collision_constant = 0.5; // avoid double counting collisions
    }
    collision_constant *= M_PI*diameter*diameter;
    collision_constant *= sqrt(8.0*zerork::KBoltzmann/(M_PI*mass));

    // the binary collision constant has units [m^3/(collisions-s)], so one
    // needs to multiply by KAvogadroNumber [# collisons/ kmol of collisons]
    // to get a constant of [m^3/(kmol-s)] - ignoring the sqrt(T) term
    collision_constant *= zerork::KAvogadroNumber;
    hard_sphere_collision_rate_constant_[step_id] = collision_constant;
  }
}


void  BinaryCollision::GetSpeciesReport(std::string *report)
{

  char format_line[zerork::MAX_FILE_LINE];

  report->clear();
  *report  = "# Binary Collision Species Report:\n";

  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   Hard sphere diameter used to estimate collision cross-sections calculated\n#   at T_ref = %.18g [K]\n",
           temperature_ref_);
  *report += std::string(format_line);
  *report += "#   Method flags: -1 == initialization\n";
  *report += "#                  1 == estimated from molecular weight\n";
  *report += "#                  2 == internal database\n";
  *report += "#                  3 == external file\n";
 
  *report += "#                                                                         [m]\n";
  *report += "#                              [K]            [m]       [Pa-s]     hard sphere\n";
  *report += "# index,  species, method, epsilon/k,        sigma,    mu(T_ref),    diameter\n"; 
  for(int j=0; j<mechp_->getNumSpecies(); ++j) {
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "%3d %14s    %3d  %9.4f  %12.4e  %12.4e  %12.4e\n",
             j,
             mechp_->getSpeciesName(j),
             method_[j],
             lennard_jones_[j].epsilon_k,
             lennard_jones_[j].sigma,
             viscosity_ref_[j],
             hard_sphere_diameter_[j]);
    *report += std::string(format_line);             
  }
}


int BinaryCollision::GetStepProbabilityReport(const double pressure,
                                              const double start_temperature,
                                              const double end_temperature,
                                              const int num_scan_points,
                                              const double min_probability,
                                              std::string *report)
{
  std::vector<double> scan_temperature;
  zerork::utilities::SortableVector initial;
  std::vector<zerork::utilities::SortableVector> step_probability;
  std::vector<double> current_step_probability;
  char format_line[zerork::MAX_FILE_LINE];
  std::string step_name;

  report->clear();
  if(num_scan_points <= 0) {
    *report  = "ERROR: In BinaryCollision::GetStepProbabilityReport(...),\n";
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "       can not generate report for %d scan temperatures\n",
             num_scan_points);
    *report += std::string(format_line);          
    return report->length();
  }
   
  scan_temperature.assign(num_scan_points,0.0);
  if(num_scan_points >= 2) {
    for(int j=0; j<num_scan_points; ++j) {
      scan_temperature[j] = start_temperature +
        (end_temperature-start_temperature)*
        static_cast<double>(j)/static_cast<double>(num_scan_points-1);
    }
  } else {
    scan_temperature[0] = start_temperature;
  }
  
  // initialize the step_probability vector<SortableVector>
  step_probability.clear();
  for(int j=0; j<num_binary_steps_; ++j) {
    initial.id = binary_step_index_[j];
    initial.v_double.assign(num_scan_points,0.0);
    step_probability.push_back(initial);
  }
  current_step_probability.assign(num_binary_steps_, 0.0);

  // compute the binary step probabilities at each scan temperature
  for(int j=0; j<num_scan_points; ++j) {
    GetStepProbability(scan_temperature[j],
                       pressure,
                       current_step_probability.data());
    for(int k=0; k<num_binary_steps_; ++k) {
      step_probability[k].v_double[j] = current_step_probability[k];
    }
  }

  // sort the binary step probabilities by the maximum probability over
  // all temperature for a give step
  SortByVectorMax(num_binary_steps_,step_probability.data());
   
   *report  = "# Binary Collision Step Probabilty Report:\n";
  if(num_scan_points >= 2) {
    // range scan
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "#   Highest probability found over T-range [%7.2f K,%7.2f K]\n",
             start_temperature,
             end_temperature);
    *report += std::string(format_line);
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "#   scanning temperature steps dT = %9.4f K\n",
             scan_temperature[1]-scan_temperature[0]);
    *report += std::string(format_line);
  } else {
    // single point 
     snprintf(format_line,zerork::MAX_FILE_LINE,
             "#   Highest probability found for T = %7.2f K\n",
	      scan_temperature[0]);   
    *report += std::string(format_line);
  }
  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   and p = %12.4e Pa (will only impact PLOG reactions)\n",
	    pressure);
  *report += std::string(format_line);

  mechp_->getReactionNameDirOfStep(step_probability[0].id, &step_name);
  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   highest probability step: %d {%s}\n",
           step_probability[0].id,
           step_name.c_str());
  *report += std::string(format_line);

  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   step %7d probability: %12.5e\n",
           step_probability[0].id,
           SortableVectorMax(&step_probability[0]));
  *report += std::string(format_line); 

  if(num_scan_points > 1) {
    WriteBinaryStepScanInBlocks(7,  // number of step probability cols per block
                                num_scan_points,
                                scan_temperature.data(),
                                &step_probability,
                                min_probability,
                                report);
  } else {
    // write out the single temperature probability as a list
    for(int j=0; j<num_binary_steps_; ++j) {

      if(step_probability[j].v_double[0] >= min_probability) {
        mechp_->getReactionNameDirOfStep(step_probability[j].id, &step_name);
        snprintf(format_line,zerork::MAX_FILE_LINE,
                 "%5d  %10.3e { %s }\n",
                 step_probability[j].id,
                 step_probability[j].v_double[0],
                 step_name.c_str());
        *report += std::string(format_line);
      } 
    } 
  }
  return report->size();
}

void BinaryCollision::WriteBinaryStepScanInBlocks(const int block_size,
                        const int num_scans,
                        const double scan_temperature[],
			const std::vector<zerork::utilities::SortableVector> *probability,
                        const double min_probability,
                        std::string *report)
{
  int step_id, list_id;
  char format_line[zerork::MAX_FILE_LINE];
  std::string step_name;
  int num_write = 0;
  int num_blocks;  // = num_binary_steps_/block_size;
  int rem_columns; // = num_binary_steps_%block_size;

  for(int j=0; j<num_binary_steps_; ++j) {
    if(SortableVectorMax(&probability->at(j)) >= min_probability) {
      ++num_write;
    } else {
      break;
    }
  }
  
  num_blocks  = num_write/block_size;
  rem_columns = num_write%block_size; 

  *report += "#   Number of bimolecur steps reported with a peak reaction\n";

  snprintf(format_line,
           zerork::MAX_FILE_LINE,
           "#   probability greater than or equal to %8.2e [-]: %d\n",
           min_probability,
           num_write);  
  (*report)+=std::string(format_line);
  
  for(int j=0; j<num_blocks; ++j) {
    // write new block header
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "# gnuplot data block (index): %d\n",
             j);
    *report += std::string(format_line);
    *report += "# column 1: [K] scan temperature\n";
    for(int k=0; k<block_size; ++k) {
      list_id = k+block_size*j;
      step_id = probability->at(list_id).id;
      mechp_->getReactionNameDirOfStep(step_id, &step_name);

      snprintf(format_line,zerork::MAX_FILE_LINE,
               "# column %d: [-] step %4d probability {max =%10.3e %s }\n",
               k+2,
               step_id,
               SortableVectorMax(&probability->at(list_id)),
               step_name.c_str());
      *report += std::string(format_line);
    }
    *report += "#------------------------------------------------------------------------------\n";
    // end new block header
    for(int k=0; k<num_scans; ++k) {
      snprintf(format_line, zerork::MAX_FILE_LINE,
               "%7.2f",
               scan_temperature[k]);
      *report += std::string(format_line);

      for(int m=0; m<block_size; m++) {
        list_id = m+block_size*j;
        snprintf(format_line,zerork::MAX_FILE_LINE,
                 " %10.3e",
                 probability->at(list_id).v_double[k]);
        *report += std::string(format_line);
      }
      *report += "\n";
    }

    if(j < num_blocks-1) {
      *report += "\n\n";
    }
  } // for loop-j over the number of complete blocks

  // complete any remaining partial blocks
  // write new block header
  if(rem_columns > 0) {
    *report += "\n\n";
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "# gnuplot data block (index): %d\n",
             num_blocks);
    *report += std::string(format_line);
    *report += "# column 1: [K] scan temperature\n";
    for(int k=0; k<rem_columns; ++k) {
      list_id = k+block_size*num_blocks;
      step_id = probability->at(list_id).id;
      mechp_->getReactionNameDirOfStep(step_id, &step_name);

      snprintf(format_line,zerork::MAX_FILE_LINE,
               "# column %d: [-] step %4d probability {max =%10.3e %s }\n",
               k+2,
               step_id,
               SortableVectorMax(&probability->at(list_id)),
               step_name.c_str());
      *report += std::string(format_line);
    }
    *report += "#------------------------------------------------------------------------------\n";
    // end new block header
    for(int k=0; k<num_scans; ++k) {
      snprintf(format_line, zerork::MAX_FILE_LINE,
               "%7.2f",
               scan_temperature[k]);
      *report += std::string(format_line);

      for(int m=0; m<rem_columns; m++) {
        list_id = m+block_size*num_blocks;
        snprintf(format_line,zerork::MAX_FILE_LINE,
                 " %10.3e",
                 probability->at(list_id).v_double[k]);
        *report += std::string(format_line);
      }
      *report += "\n";
    }
  } // end if(rem_columns > 0)
}

} // end of namespace transport
