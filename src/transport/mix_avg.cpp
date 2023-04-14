#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <vector> // needed for local string utilities
#include <string> // needed for local string utilities

#include <utilities/string_utilities.h>
#include <utilities/file_utilities.h> // also brings in string utilities

#include "mix_avg.h"


namespace transport
{

MixAvg::MixAvg()
{
  initialized_ = false;
  num_species_ = -1;
  mechanism_ = NULL;
  mechanism_owner_ = true; 
}

MixAvg::~MixAvg()
{
  if(mechanism_owner_ && mechanism_ != NULL) {
    delete mechanism_;
  }
}

int MixAvg::Initialize(const std::vector<std::string> &input_files,
                       const std::string &log_name,
                       const double conductivity_multiplier)
{
  int flag;
  FILE *log_fptr = fopen(log_name.c_str(),"a");

  if(log_fptr == NULL) {
    printf("# ERROR: In MixAvg::Initialize(),\n"
           "#        log file named = %s could not be opened.\n",
           log_name.c_str());
    fflush(stdout);
    return LOG_FILE_ERROR;
  }
  log_name_ = std::string(log_name);

  if(input_files.size() != 3) {
    fprintf(log_fptr,
            "# ERROR: In MixAvg::Initialize(),\n"
            "#        input_files vector incorrect size=%d.\n"
            "#        Correct size is three with the following definition:\n"
            "#            input_files[0] = mechanism file\n"
            "#            input_files[1] = thermodynamics file\n"
            "#            input_files[2] = transport file\n",
            (int)input_files.size());

    fflush(log_fptr);
    fclose(log_fptr);
    return INPUT_FILE_ERROR;
  }
  fclose(log_fptr);

  mechanism_ = new zerork::mechanism(input_files[0].c_str(),
                                     input_files[1].c_str(),
                                     log_name_.c_str());

  if(mechanism_ == NULL) {
    log_fptr = fopen(log_name.c_str(),"a");
    fprintf(log_fptr,
            "# ERROR: In MixAvg::Initialize(),\n"
            "#        could not instantiate zerork::mechanism from\n"
            "#            mechanism      file = %s\n"
            "#            thermodynamics file = %s\n",
            input_files[0].c_str(),
            input_files[1].c_str());
    fflush(log_fptr);
    fclose(log_fptr);
    return INPUT_FILE_ERROR;
  }

  // set up the molecular mass arrays
  num_species_ = mechanism_->getNumSpecies();
  molecular_mass_.assign(num_species_, 0.0);
  inv_molecular_mass_.assign(num_species_, 0.0);
  mechanism_->getMolWtSpc(&molecular_mass_[0]);
  for(int j=0; j<num_species_; ++j) {
    inv_molecular_mass_[j] = 1.0/molecular_mass_[j];
  }
  species_workspace_.assign(num_species_, 0.0);

  lewisnumbers.assign(num_species_, 0.0);
  Dmass.assign(num_species_, 0.0);
  invDij.assign(num_species_*num_species_, 0.0);

  multiplier_ = conductivity_multiplier;

  // parse the transport file for species transport parameters
  std::string error_message;
  flag = ParseTransportFile(input_files[2],
                            "!", // comment character(s)
                            &error_message);
  if(flag == NO_ERROR) {
    initialized_ = true;
  } else {
    initialized_ = false;
  }

  return flag;
}

int MixAvg::Initialize(zerork::mechanism* mechanism,
                       const std::vector<std::string> &input_files,
                       const std::string &log_name,
                       const double conductivity_multiplier)
{
  int flag;
  FILE *log_fptr = fopen(log_name.c_str(),"a");

  if(log_fptr == NULL) {
    printf("# ERROR: In MixAvg::Initialize(),\n"
           "#        log file named = %s could not be opened.\n",
           log_name.c_str());
    fflush(stdout);
    return LOG_FILE_ERROR;
  }
  log_name_ = std::string(log_name);

  if(input_files.size() != 1) {
    fprintf(log_fptr,
            "# ERROR: In MixAvg::Initialize(),\n"
            "#        input_files vector incorrect size=%d.\n"
            "#        Correct size is one with the following definition:\n"
            "#            input_files[0] = transport file\n",
            (int)input_files.size());

    fflush(log_fptr);
    fclose(log_fptr);
    return INPUT_FILE_ERROR;
  }
  fclose(log_fptr);

  mechanism_ = mechanism;
  mechanism_owner_ = false;

  if(mechanism_ == NULL) {
    log_fptr = fopen(log_name.c_str(),"a");
    fprintf(log_fptr,
            "# ERROR: In MixAvg::Initialize(),\n"
            "#        could not instantiate zerork::mechanism from\n"
            "#            mechanism      file = %s\n"
            "#            thermodynamics file = %s\n",
            input_files[0].c_str(),
            input_files[1].c_str());
    fflush(log_fptr);
    fclose(log_fptr);
    return INPUT_FILE_ERROR;
  }

  // set up the molecular mass arrays
  num_species_ = mechanism_->getNumSpecies();
  molecular_mass_.assign(num_species_, 0.0);
  inv_molecular_mass_.assign(num_species_, 0.0);
  mechanism_->getMolWtSpc(&molecular_mass_[0]);
  for(int j=0; j<num_species_; ++j) {
    inv_molecular_mass_[j] = 1.0/molecular_mass_[j];
  }
  species_workspace_.assign(num_species_, 0.0);

  lewisnumbers.assign(num_species_, 0.0);
  Dmass.assign(num_species_, 0.0);
  invDij.assign(num_species_*num_species_, 0.0);

  multiplier_ = conductivity_multiplier;

  // parse the transport file for species transport parameters
  std::string error_message;
  flag = ParseTransportFile(input_files[0],
                            "!", // comment character(s)
                            &error_message);
  if(flag == NO_ERROR) {
    initialized_ = true;
  } else {
    initialized_ = false;
  }

  return flag;
}

int MixAvg::GetMixtureViscosity(const MassTransportInput &input,
                                       double *viscosity) const
{
  int flag = GetSpeciesViscosity(input, &species_workspace_[0]);

  if(flag == NO_ERROR) {

    const int num_species = num_species_;
    double molecular_mass_mix;
    double inv_molecular_mass_mix = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;

    for(int j=0; j<num_species; ++j) {

      double product = input.mass_fraction_[j]*inv_molecular_mass_[j];

      inv_molecular_mass_mix += product;
      sum1 += product*species_workspace_[j];
      // TODO: Add vectorized 1/x operation if major cost
      sum2 += product/species_workspace_[j];
    }
    molecular_mass_mix = 1.0/inv_molecular_mass_mix;
    *viscosity = 0.5*(molecular_mass_mix*sum1 +
                      1.0/(molecular_mass_mix*sum2));

  }
  return flag;
}
int MixAvg::GetSpeciesViscosity(const MassTransportInput &input,
                                       double *viscosity) const
{
  if(input.temperature_ > 0.0) {

    // TODO: add vectorization if needed
    const int num_species = num_species_;
    const double temperature = input.temperature_;

    for(int j=0; j<num_species; ++j) {
      viscosity[j] = mucoeff_[j]*sqrt(temperature)/omega_mu(temperature*kOverEps_[j]);
    }

    return NO_ERROR;
  }

  return NON_POSITIVE_TEMPERATURE;
}

int MixAvg::GetMixtureConductivity(const MassTransportInput &input,
                                          double *conductivity) const
{
  int flag = GetSpeciesConductivity(input, &species_workspace_[0]);

  if(flag == NO_ERROR) {

    const int num_species = num_species_;
    double molecular_mass_mix;
    double inv_molecular_mass_mix = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;

    for(int j=0; j<num_species; ++j) {

      double product = input.mass_fraction_[j]*inv_molecular_mass_[j];

      inv_molecular_mass_mix += product;
      sum1 += product*species_workspace_[j];
      // TODO: Add vectorized 1/x operation if major cost
      sum2 += product/species_workspace_[j];
    }
    molecular_mass_mix = 1.0/inv_molecular_mass_mix;
    *conductivity = 0.5*(molecular_mass_mix*sum1 +
                         1.0/(molecular_mass_mix*sum2))*multiplier_;

  }
  return flag;
}

int MixAvg::GetSpeciesConductivity(const MassTransportInput &input,
                                          double *conductivity) const
{
  if(input.temperature_ > 0.0) {

    // TODO: add vectorization if needed
    const int num_species = num_species_;
    const double temperature = input.temperature_;
    const double sqrt_temperature = sqrt(temperature);
    double beta;
    double omegamu;
    std::vector<double> Cp_sp(num_species);
    const double gasConstant = 8314.462;

    const double specific_heat_cp_mass =
      mechanism_->getMassCpFromTY(input.temperature_,
                                  &input.mass_fraction_[0],
                                  &Cp_sp[0]);


    for(int j=0; j<num_species; ++j) {
      // Modified Eucken formula
      omegamu = omega_mu(temperature*kOverEps_[j]);
      beta = 1.2*omegamu/omega_D(temperature*kOverEps_[j]);
      conductivity[j] = mucoeff_[j]*sqrt_temperature/omegamu*
	(beta*Cp_sp[j] + (3.75-2.5*beta)*gasConstant*inv_molecular_mass_[j]);
    }

    return NO_ERROR;
  }

  return NON_POSITIVE_TEMPERATURE;
}

int MixAvg::GetSpeciesMassFlux(const MassTransportInput &input,
			       const size_t ld_species_mass_flux,
			       double *conductivity_mix,
			       double *specific_heat_mix,
			       double *species_mass_flux,
			       double *species_lewis_numbers) const
{
  int flag = NO_ERROR;
  double conductivity_mix_local=0.0;
  if(conductivity_mix != nullptr && *conductivity_mix > 0) {
      conductivity_mix_local = *conductivity_mix;
  } else {
    flag = GetMixtureConductivity(input,
                                  &conductivity_mix_local);
    if(conductivity_mix != nullptr) {
      *conductivity_mix = conductivity_mix_local;
    }
  }
  const double gasConstant = 8314.46;

  if(flag == NO_ERROR) {

    const int num_species = num_species_;
    const size_t num_dimensions = input.num_dimensions_;

    const double molecular_mass_mix =
      mechanism_->getMolWtMixFromY(&input.mass_fraction_[0]);
    const double inv_molecular_mass_mix = 1.0/molecular_mass_mix;

    double specific_heat_cp_mass = 0;
    if(specific_heat_mix != nullptr && *specific_heat_mix > 0) {
        specific_heat_cp_mass = *specific_heat_mix;
    } else {
      specific_heat_cp_mass = mechanism_->getMassCpFromTY(input.temperature_,
                                                          &input.mass_fraction_[0]);
      if(specific_heat_mix != nullptr) {
        *specific_heat_mix = specific_heat_cp_mass;
      }
    }

    // zero the species mass flux
    size_t flux_id = 0;
    size_t grad_id = 0;
    for(size_t j=0; j<num_dimensions; ++j) {

      for(int k=0; k<num_species; ++k) {
        species_mass_flux[flux_id] = 0.0;
        ++flux_id;
      }
      flux_id += ld_species_mass_flux;
    }

    std::vector<double> mass_flux_sum;
    mass_flux_sum.assign(num_dimensions,0.0);

    // compute molecular_mass_mix*\sum_i (1/molecular_mass[i])*
    //                                    \grad(mass_fraction[i])
    grad_id = 0;
    for(size_t j=0; j<num_dimensions; ++j) {

      for(int k=0; k<num_species; ++k) {
        mass_flux_sum[j] +=
          input.grad_mass_fraction_[grad_id]*inv_molecular_mass_[k];
	++grad_id;
      } //k species loop
      grad_id += input.ld_grad_mass_fraction_;
      mass_flux_sum[j] *= molecular_mass_mix;
    } //j dimensions loop

    flux_id = 0;
    grad_id = 0;

    for(size_t j=0; j<num_dimensions; ++j) {

      const double dcoeff = 419.75742*input.pressure_/ sqrt(pow(input.temperature_,3.0)*1000);
      const double rho = input.pressure_*molecular_mass_mix/(gasConstant*input.temperature_);

      // Compute Mass diffusion term Dmass_k
      for(int k=0; k<num_species; ++k) {

	double num = 0.0;
	double den = 0.0;
	for(int l=0; l<num_species; ++l) {

	  invDij[k*num_species + l] = dcoeff*diam2_[k*num_species+l]*
	    omega_D(input.temperature_*sqrtkOverEps_[k]*sqrtkOverEps_[l])*
	    sqrtmass_[k*num_species+l];

	  if(l != k) {
	    num += input.mass_fraction_[l];
	    den += input.mass_fraction_[l]*inv_molecular_mass_[l]*invDij[k*num_species + l];
	  }

	} // l species loop
        if(den != 0.0) {
          Dmass[k] = rho*num*inv_molecular_mass_mix/den;
        } else {
          Dmass[k] = rho/invDij[k*num_species + k];
        }
	// Compute Lewis numbers from DMass_k
	lewisnumbers[k] = conductivity_mix_local/(specific_heat_cp_mass*Dmass[k]);

      } //k species loop

      // compute the uncorrected species mass flux
      for(int k=0; k<num_species; ++k) {

	species_mass_flux[flux_id] =
          input.grad_mass_fraction_[grad_id] -
          input.mass_fraction_[k]*mass_flux_sum[j];

	// Compute species flux  -- mass
	species_mass_flux[flux_id] *= -Dmass[k];
	// Save Lewis numbers
	species_lewis_numbers[flux_id] = lewisnumbers[k];

        ++flux_id;
        ++grad_id;
      } //species loop
      flux_id += ld_species_mass_flux;
      grad_id += input.ld_grad_mass_fraction_;

    } //dimension loop

    // compute the species mass flux sum for correction
    mass_flux_sum.assign(num_dimensions,0.0);
    flux_id = 0;

    for(size_t j=0; j<num_dimensions; ++j) {

      for(int k=0; k<num_species; ++k) {
	mass_flux_sum[j] += species_mass_flux[flux_id];
        ++flux_id;
      }
      flux_id += ld_species_mass_flux;
    }

    // apply the species mass flux sum correction
    flux_id = 0;

    for(size_t j=0; j<num_dimensions; ++j) {

      double correction = mass_flux_sum[j];

      for(int k=0; k<num_species; ++k) {

	species_mass_flux[flux_id] -= input.mass_fraction_[k]*correction;
        ++flux_id;
      }
      flux_id += ld_species_mass_flux;
    }

  }

  return flag;
}


int MixAvg::GetSpeciesMassFluxFrozenThermo(const MassTransportInput &input,
					   const size_t ld_species_mass_flux,
					   double *conductivity_mix,
					   double *specific_heat_mix,
					   double *species_mass_flux,
					   double *species_lewis_numbers) const
{
  int flag = NO_ERROR;
  double conductivity_mix_local=0.0;
  if(conductivity_mix != nullptr && *conductivity_mix > 0) {
      conductivity_mix_local = *conductivity_mix;
  } else {
    flag = GetMixtureConductivity(input,
                                  &conductivity_mix_local);
    if(conductivity_mix != nullptr) {
      *conductivity_mix = conductivity_mix_local;
    }
  }
  double dcoeff = 0.0;
  double rho = 0.0;
  const double gasConstant = 8314.46;

  if(flag == NO_ERROR) {

    const int num_species = num_species_;
    const size_t num_dimensions = input.num_dimensions_;

    const double molecular_mass_mix =
      mechanism_->getMolWtMixFromY(&input.mass_fraction_[0]);
    const double inv_molecular_mass_mix = 1.0/molecular_mass_mix;

    double specific_heat_cp_mass = 0;
    if(specific_heat_mix != nullptr && *specific_heat_mix > 0) {
        specific_heat_cp_mass = *specific_heat_mix;
    } else {
      specific_heat_cp_mass = mechanism_->getMassCpFromTY(input.temperature_,
                                                          &input.mass_fraction_[0]);
      if(specific_heat_mix != nullptr) {
        *specific_heat_mix = specific_heat_cp_mass;
      }
    }

    // zero the species mass flux
    size_t flux_id = 0;
    size_t grad_id = 0;
    for(size_t j=0; j<num_dimensions; ++j) {

      for(int k=0; k<num_species; ++k) {
        species_mass_flux[flux_id] = 0.0;
        ++flux_id;
      }
      flux_id += ld_species_mass_flux;
    }

    std::vector<double> mass_flux_sum;
    mass_flux_sum.assign(num_dimensions,0.0);

    // compute molecular_mass_mix*\sum_i (1/molecular_mass[i])*
    //                                    \grad(mass_fraction[i])
    grad_id = 0;
    for(size_t j=0; j<num_dimensions; ++j) {

      for(int k=0; k<num_species; ++k) {
        mass_flux_sum[j] +=
          input.grad_mass_fraction_[grad_id]*inv_molecular_mass_[k];
	++grad_id;
      } //k species loop
      grad_id += input.ld_grad_mass_fraction_;
      mass_flux_sum[j] *= molecular_mass_mix;
    } //j dimensions loop

    flux_id = 0;
    grad_id = 0;

    for(size_t j=0; j<num_dimensions; ++j) {

      dcoeff = 419.75742*input.pressure_/ sqrt(pow(input.temperature_,3.0)*1000);
      rho = input.pressure_*molecular_mass_mix/(gasConstant*input.temperature_);

      // Compute Mass diffusion term Dmass_k
      for(int k=0; k<num_species; ++k) {

	double num = 0.0;
	double den = 0.0;
	for(int l=0; l<num_species; ++l) {

	  invDij[k*num_species + l] = dcoeff*diam2_[k*num_species+l]*
	    omega_D(input.temperature_*sqrtkOverEps_[k]*sqrtkOverEps_[l])*
	    sqrtmass_[k*num_species+l];

	  if(l != k) {
	    num += input.mass_fraction_[l];
	    den += input.mass_fraction_[l]*inv_molecular_mass_[l]*invDij[k*num_species + l];
	  }

	} // l species loop
        if(den != 0.0) {
          Dmass[k] = rho*num*inv_molecular_mass_mix/den;
        } else {
          Dmass[k] = rho/invDij[k*num_species + k];
        }
	// Compute Lewis numbers from DMass_k
	lewisnumbers[k] = conductivity_mix_local/(specific_heat_cp_mass*Dmass[k]);

      } //k species loop

      // compute the uncorrected species mass flux
      for(int k=0; k<num_species; ++k) {

	species_mass_flux[flux_id] =
          input.grad_mass_fraction_[grad_id] -
          input.mass_fraction_[k]*mass_flux_sum[j];

	// Compute species flux  -- mass
	species_mass_flux[flux_id] *= -Dmass[k];
	// Save Lewis numbers
	species_lewis_numbers[flux_id] = lewisnumbers[k];

        ++flux_id;
        ++grad_id;
      } //species loop
      flux_id += ld_species_mass_flux;
      grad_id += input.ld_grad_mass_fraction_;

    } //dimension loop

  }

  return flag;
}


// The format of the transport file
int MixAvg::ParseTransportFile(const std::string &transport_file,
                                      const std::string &comment_chars,
                                      std::string *error_message)
{
  const int num_species = mechanism_->getNumSpecies();
  bool incomplete_file = false;
  int line_num = 0;
  size_t comment_pos;
  size_t num_fields;
  std::string line, sub_line;
  std::vector<std::string> fields;
  std::vector<int> species_line_num;
  std::map<std::string, int> species_id_from_name;
  std::map<std::string, int>::const_iterator map_iter;
  std::ifstream input_file(transport_file.c_str());

  FILE *log_fptr = fopen(log_name_.c_str(),"a");

  error_message->clear();

  if(input_file) {

    // build map from species names in mechanism
    for(int j=0; j<num_species; ++j) {
      std::string species_name = std::string(mechanism_->getSpeciesName(j));
      species_id_from_name[species_name] = j;
    }

    // initialize data parameters
    shape_.assign(num_species, 0);
    kOverEps_.assign(num_species, 0.0);
    sqrtkOverEps_.assign(num_species, 0.0);
    sigma_.assign(num_species, 0.0);
    mu_.assign(num_species, 0.0);
    alpha_.assign(num_species, 0.0);
    Zrot_.assign(num_species, 0.0);

    sqrtmass_.assign(num_species*num_species, 0.0);
    diam2_.assign(num_species*num_species, 0.0);

    mucoeff_.assign(num_species, 0.0);
    species_line_num.assign(num_species, -1);

    while(zerork::utilities::GetAnyLine(input_file, &line)) {

      ++line_num;

      // extract the portion of the line before any comment
      comment_pos = line.find_first_of(comment_chars);
      sub_line.assign(line, 0, comment_pos);
      num_fields = zerork::utilities::SplitStringToVector(sub_line,
                                                  zerork::utilities::WHITESPACE,
                                                  &fields);
      if(num_fields >= 7) {
        map_iter = species_id_from_name.find(fields[0]);

        if(map_iter != species_id_from_name.end()) {
          // found species name
          int species_id = map_iter->second;

          // only record data from the first appearance, otherwise ignore
          // a repeated entry
          if(species_line_num[species_id] <= 0) {
            species_line_num[species_id] = line_num;

            shape_[species_id]        = atof(fields[1].c_str());
            kOverEps_[species_id]     = 1.0/atof(fields[2].c_str());
            sigma_[species_id]        = atof(fields[3].c_str());
            mu_[species_id]           = atof(fields[4].c_str());
            alpha_[species_id]        = atof(fields[5].c_str());
            Zrot_[species_id]         = atof(fields[6].c_str());

	    sqrtkOverEps_[species_id] = sqrt(kOverEps_[species_id]);
	    mucoeff_[species_id] = 2.6693e-6 * sqrt(molecular_mass_[species_id])/pow(sigma_[species_id],2.0);
	  }
          else if(log_fptr != NULL) {
            fprintf(log_fptr,
                    "# INFO: At line %d in transport file %s,\n"
                    "#       skipping duplicate data for species %s that was\n"
                    "#       previously specified at line %d.\n",
                    line_num,
                    transport_file.c_str(),
                    fields[0].c_str(),
                    species_line_num[species_id]);
            fflush(log_fptr);

          }
        } // end if(map_iter != species_id_from_name.end())
        if(map_iter == species_id_from_name.end() &&
           log_fptr != NULL) {

          fprintf(log_fptr,
                  "# INFO: At line %d in transport file %s,\n"
                  "#       skipping data for species %s that is not in the\n"
                  "#       mechanism file.\n",
                  line_num,
                  transport_file.c_str(),
                  fields[0].c_str());
          fflush(log_fptr);

        }
      } // end if(num_fields >= 7)
    } // end while(getline(input_file,line))

    // Check that the transport parameters were defined for all species by
    // checking that the line numbers all greater than zero
    incomplete_file = false;
    for(int j=0; j<num_species; ++j) {
      if(species_line_num[j] <= 0) {
        incomplete_file = true;
        if(error_message->size() > 0) {
          error_message->append("\n");
        }
        error_message->append("No transport parameters defined for species " +
                              std::string(mechanism_->getSpeciesName(j)) +
                              " in transport file = " +
                              transport_file + ".");
      }
    }
    if(incomplete_file) {

      if(log_fptr != NULL) {
        fprintf(log_fptr,"%s\n",error_message->c_str());
        fclose(log_fptr);
      }
      return INPUT_FILE_ERROR;
    }

    if(log_fptr != NULL) {
      fclose(log_fptr);
    }

    for(int k=0; k<num_species; ++k) {
      for(int l=0; l<num_species; ++l) {
	sqrtmass_[k*num_species+l] = sqrt(molecular_mass_[k]*molecular_mass_[l] /
					 (molecular_mass_[k] + molecular_mass_[l]) );
	diam2_[k*num_species+l] = pow(sigma_[k] + sigma_[l],2.0);
      }
    }


    // TODO: add transport to log file
    return NO_ERROR;

  } else {

    error_message->assign("Can not open transport file = " +
                          transport_file + ".");
    if(log_fptr != NULL) {
      fprintf(log_fptr,"%s\n",error_message->c_str());
      fclose(log_fptr);
    }
    return INPUT_FILE_ERROR;
  } // if input file
}

double MixAvg::omega_D(double t) const
{
  static double m1 = 6.8728271691;
  static double m2 = 9.4122316321;
  static double m3 = 7.7442359037;
  static double m4 = 0.23424661229;
  static double m5 = 1.45337701568;
  static double m6 = 5.2269794238;
  static double m7 = 9.7108519575;
  static double m8 = 0.46539437353;
  static double m9 = 0.00041908394781;

  double num = m1 + t * (m2 + t * (m3 + t * m4));
  double den = m5 + t * (m6 + t * (m7 + t * (m8 + t * m9)));
  return num / den;
}

  double MixAvg::omega_mu (double t) const
  {
    static double m1 = 3.3530622607;
    static double m2 = 2.53272006;
    static double m3 = 2.9024238575;
    static double m4 = 0.11186138893;
    static double m5 = 0.8662326188;
    static double m6 = 1.3913958626;
    static double m7 = 3.158490576;
    static double m8 = 0.18973411754;
    static double m9 = 0.00018682962894;

    double num = m1 + t * (m2 + t * (m3 + t * m4));
    double den = m5 + t * (m6 + t * (m7 + t * (m8 + t * m9)));
    return num / den;
  }


} // namespace transport
