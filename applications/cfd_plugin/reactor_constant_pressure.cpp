#include <algorithm> //std::max,min

#include "reactor_constant_pressure.h"

#include "nvector/nvector_serial.h"

ReactorConstantPressure::ReactorConstantPressure(std::shared_ptr<zerork::mechanism> mech_ptr)
 :
    ReactorNVectorSerial(mech_ptr)
{}

ReactorConstantPressure::~ReactorConstantPressure()
{}

void ReactorConstantPressure::InitializeState(
    const double reactor_time,
    const int n_reactors,
    const double *T,
    const double *P,
    const double *mf,
    const double *dpdt,
    const double *e_src,
    const double *y_src)
{
  assert(n_reactors == 1);
  ReactorNVectorSerial::Reset(); //re-initialize to get deterministic behavior across solves
  pressure_ = *P;
  initial_temperature_ = *T;
  inverse_density_ = 1.0/mech_ptr_->getDensityFromTPY(*T,*P,mf);
  initial_energy_ = mech_ptr_->getMassEnthalpyFromTY(initial_temperature_,mf);
  double *y_ptr = NV_DATA_S(state_);
  for(int k = 0; k < num_species_; ++k) {
    y_ptr[k] = mf[k];
  }
  y_ptr[num_species_] = *T/double_options_["reference_temperature"];
  dpdt_ = *dpdt;
  e_src_ = *e_src;
  y_src_ = y_src;
}

void ReactorConstantPressure::GetState(
    double *T,
    double *P,
    double *mf)
{
  *T = NV_Ith_S(state_,num_species_)*double_options_["reference_temperature"];
  double *y_ptr = NV_DATA_S(state_);
  for(int k = 0; k < num_species_; ++k) {
    mf[k] = y_ptr[k];
  }
  *P = pressure_;
}

int ReactorConstantPressure::GetTimeDerivative(const double reactor_time,
                                                N_Vector state,
                                                N_Vector derivative)
{
  //double startTime=getHighResolutionTime();
  double * y_ptr = NV_DATA_S(state);
  double * ydot_ptr = NV_DATA_S(derivative);
  double * net_production_rates_ptr = &(net_production_rates_[0]);
  double * creation_rates_ptr = &(creation_rates_[0]);
  double * destruction_rates_ptr = &(destruction_rates_[0]);
  double * forward_rates_of_production_ptr = &(forward_rates_of_production_[0]);
  double * energy_ptr = &(energy_[0]);
  double * cp_mass_ptr = &(cx_mass_[0]);
  const int num_spec = num_species_;
  double temperature = y_ptr[num_spec]*double_options_["reference_temperature"];
  if(temperature <= 0.0) return 1;
#define TLIMIT 1.0e4
  temperature = std::min(temperature,TLIMIT);
  //TODO:
  //check_negative_fractions(y_ptr)
  //temperature = clip_temperature(temperature);

//  static int ncalls = 0;
//  ncalls++;
//  if(ncalls > 20) { exit(1); }
//  printf("%d - %g, %g\n",ncalls, reactor_time, temperature);

//#define ATOL_SAFETY_FACTOR 1.0e10
//  const double check_val = -double_options_["abs_tol"]*ATOL_SAFETY_FACTOR;
//  for(int j=0; j<num_species_; ++j) {
//    if(y_ptr[j] < check_val) {
//      //if(int_options_["verbosity"] > 4) { 
//        printf("WARNING: Significantly negative species in "
//               "Zero-RK RHS.  Consider reducing tolerances.\n");
//      //}
//      return 1;
//    }
//  }
  if(e_src_ != 0) {
    double energy = initial_energy_ + e_src_*reactor_time; //N.B. assumes t_0 = 0.0
    temperature = mech_ptr_->getTemperatureFromHY(energy, y_ptr, temperature);
  }
  mech_ptr_->getEnthalpy_RT(temperature,energy_ptr);

  inverse_density_ = 1.0/mech_ptr_->getDensityFromTPY(temperature,pressure_,y_ptr);

  mean_cx_mass_ = mech_ptr_->getMassCpFromTY(temperature,y_ptr,cp_mass_ptr);

  // set concentration via density and mass fraction
  mech_ptr_->getCfromVY(inverse_density_,y_ptr,&concentrations_[0]);

  // compute the molar production rates at the current state_ (aka wdot)
  mech_ptr_->getReactionRates(temperature, &concentrations_[0],
                         net_production_rates_ptr, creation_rates_ptr, destruction_rates_ptr,
                         forward_rates_of_production_ptr);

  // ydot = [kmol/m^3/s] * [kg/kmol] * [m^3/kg] = [(kg spec j)/(kg mix)/s]
  if(y_src_ != nullptr) {
    for(int j=0; j<num_spec; ++j) {
      ydot_ptr[j]=(net_production_rates_[j]*mol_wt_[j])*inverse_density_ + y_src_[j];
    }
  } else {
    for(int j=0; j<num_spec; ++j) {
      ydot_ptr[j]=(net_production_rates_[j]*mol_wt_[j])*inverse_density_;
    }
  }

  double energy_sum=0.0;
  for(int j=0; j<num_spec; ++j) {
    energy_sum += energy_ptr[j]*net_production_rates_ptr[j];
  }
  double dT_dt = -energy_sum * mech_ptr_->getGasConstant() *
                 y_ptr[num_spec] * inverse_density_ / mean_cx_mass_;
  ydot_ptr[num_spec]= dT_dt + e_src_ / mean_cx_mass_;

  //num_function_calls += 1;
  //function_time += getHighResolutionTime() - startTime;
  return 0;
}

