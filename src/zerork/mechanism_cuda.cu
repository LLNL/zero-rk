#include "mechanism_cuda.h"
#include "zerork_cuda_defs.h"
#include "mechanism_kernels.h"
#include "nasa_poly_cuda.h"
#include "rate_const_cuda.h"
#include "perf_net_cuda.h"
#include <cassert>

#include <thrust/device_vector.h>


namespace zerork {

mechanism_cuda::mechanism_cuda(const char *mechFileName,
			     const char *thermFileName, 
			     const char *convertFileName,
                             int verbosity_inp,
                             int nReactorsMax)
  : mechanism(mechFileName,thermFileName,convertFileName, verbosity_inp)
{
  assert(nReactorsMax>0);
  allocate_device_mem();
  delete thermo;
  delete infoNet;
  delete Kconst;
  delete perfNet;

  ckr::CKReader ckrobj;
 
  //We alread parsed once successfully in mechanism.
  //  Don't output again and don't bother with error messages
  //  so we don't have to include the logic for when to print
  //  error messages.
  std::string ckFileStr = "";

  bool parse_ok = ckrobj.read(mechFileStr,thermFileStr,ckFileStr);
  assert(parse_ok);//already parsed in non cuda mech

  ckrobj.verbose=false;
  for(int j=0; j<nSpc; j++)
  {
      // set Species.index to its position in the species list and the species
      // name map
      ckrobj.species[j].index=j;
      ckrobj.speciesData[ckrobj.species[j].name].index=j;
  }

  initialize_ptrs_cuda(&ckrobj,nReactorsMax);
}

mechanism_cuda::~mechanism_cuda()
{
  free_device_mem();
}


void mechanism_cuda::allocate_device_mem()
{
  //mechanism constructed, so performance
  //arrays should be calculated.
  checkCudaError
  (
      cudaMalloc((void**)&molWt_dev,sizeof(double)*nSpc),
      "cudaMalloc(... molWt_dev ...)"
  );
  cudaMemcpy(molWt_dev,molWt,sizeof(double)*nSpc,cudaMemcpyHostToDevice);
  checkCudaError
  (
      cudaMalloc((void**)&invMolWt_dev,sizeof(double)*nSpc),
      "cudaMalloc(... invMolWt_dev ...)"
  );
  cudaMemcpy(invMolWt_dev,invMolWt,sizeof(double)*nSpc,cudaMemcpyHostToDevice);
  checkCudaError
  (
      cudaMalloc((void**)&RuInvMolWt_dev,sizeof(double)*nSpc),
      "cudaMalloc(... RuInvMolWt_dev ...)"
  );
  cudaMemcpy(RuInvMolWt_dev,RuInvMolWt,sizeof(double)*nSpc,cudaMemcpyHostToDevice);
}


void mechanism_cuda::free_device_mem()
{
  cudaFree(molWt_dev);
  cudaFree(invMolWt_dev);
  cudaFree(RuInvMolWt_dev);
}

//Same as mechanism for now
void mechanism_cuda::initialize_ptrs_cuda(ckr::CKReader *ckrobj, int nReactorsMax)
{
  int j,k;
  int nLoCoef,nHiCoef;
  double *Tlo, *Thi, *coef;

  // fill and create the data for the thermodynamics object
  coef = new double[nSpc*NUM_THERMO_POLY_D5R2];
  Tlo  = new double[nSpc];
  Thi  = new double[nSpc];
  
  for(j=0; j<nSpc; j++)
  {
    Tlo[j]=ckrobj->species[j].tlow;
    Thi[j]=ckrobj->species[j].thigh;
    coef[j*NUM_THERMO_POLY_D5R2+0]=ckrobj->species[j].tmid;

    nLoCoef=ckrobj->species[j].lowCoeffs.size();
//    if(nLoCoef != NUM_COEF_RANGE)
//    {
//      printf("ERROR: unexpected number of low temperature coeffients\n");
//      printf("       for species %d (%s) - %d coeffiecients in CKReader\n",
//             j+1,ckrobj->species[j].name.c_str(),nLoCoef);
//    }
    assert(nLoCoef == NUM_COEF_RANGE);
    for(k=0; k<NUM_COEF_RANGE; k++)
    {
      coef[j*NUM_THERMO_POLY_D5R2+k+1]= ckrobj->species[j].lowCoeffs[k];
    }
    nHiCoef=ckrobj->species[j].highCoeffs.size();
//    if(nHiCoef != NUM_COEF_RANGE)
//    {
//      printf("ERROR: unexpected number of high temperature coeffients\n");
//      printf("       for species %d (%s) - %d coeffiecients in CKReader\n",
//             j+1,ckrobj->species[j].name.c_str(),nHiCoef);
//    }
    assert(nHiCoef == NUM_COEF_RANGE);
    for(k=0; k<NUM_COEF_RANGE; k++)
    {
      coef[j*NUM_THERMO_POLY_D5R2+k+1+NUM_COEF_RANGE]=
          ckrobj->species[j].highCoeffs[k];
    }
  }

  thermo = new nasa_poly_group_cuda(nSpc,coef,Tlo,Thi);
  infoNet = new info_net(ckrobj);
  Kconst = new rate_const_cuda(ckrobj,infoNet,thermo,nReactorsMax);
  perfNet = new perf_net_cuda(*infoNet,*Kconst);

  delete [] Thi;
  delete [] Tlo;
  delete [] coef;
}

int mechanism_cuda::nReactorsMax()
{
  return static_cast<rate_const_cuda*>(Kconst)->nReactorsMax();
}


void mechanism_cuda::getReactionRates_CUDA_mr_dev(const int nReactors, const double T_dev[], const double C_dev[],
                                     double netOut_dev[], double createOut_dev[],
                                     double destroyOut_dev[], double stepOut_dev[])
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  static_cast<perf_net_cuda*>(perfNet)->calcRatesFromTC_CUDA_mr_dev(nReactors,T_dev,C_dev, nullptr, netOut_dev,createOut_dev,destroyOut_dev,stepOut_dev);
}

void mechanism_cuda::getReactionRatesLimiter_CUDA_mr_dev(const int nReactors, const double T_dev[], const double C_dev[], const double stepLimiter_dev[],
                                     double netOut_dev[], double createOut_dev[],
                                     double destroyOut_dev[], double stepOut_dev[])
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  static_cast<perf_net_cuda*>(perfNet)->calcRatesFromTC_CUDA_mr_dev(nReactors,T_dev,C_dev,stepLimiter_dev,netOut_dev,createOut_dev,destroyOut_dev,stepOut_dev);
}


void mechanism_cuda::getMassCpFromTY_mr_dev(const int nReactors, const double T_dev[], const double y_dev[],
                                      double cpSpc_dev[], double cpReactors_dev[]) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  cudaMemsetAsync(cpReactors_dev,0,sizeof(double)*nReactors, (cudaStream_t) 0);
  static_cast<nasa_poly_group_cuda*>(thermo)->getCp_R_CUDA_mr(nReactors,T_dev,cpSpc_dev, (cudaStream_t) 0);
  //Convert to mass based Cp and calculate means
  meanCpMass_mr(nReactors,nSpc,RuInvMolWt_dev,y_dev,cpSpc_dev,cpReactors_dev);
}


void mechanism_cuda::getMassCvFromTY_mr_dev(const int nReactors, const double T_dev[], const double y_dev[],
                                      double cvSpc_dev[], double cvReactors_dev[]) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  cudaMemsetAsync(cvReactors_dev,0,sizeof(double)*nReactors, (cudaStream_t) 0);
  static_cast<nasa_poly_group_cuda*>(thermo)->getCp_R_CUDA_mr(nReactors,T_dev,cvSpc_dev, (cudaStream_t) 0);
  Cv_from_Cp_mr(nReactors,nSpc,RuInvMolWt_dev,y_dev,cvSpc_dev,cvReactors_dev);
}


void mechanism_cuda::getEnthalpy_RT_mr_dev(const int nReactors, const double T_dev[], double h_RT_dev[]) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  static_cast<nasa_poly_group_cuda*>(thermo)->getH_RT_CUDA_mr(nReactors,T_dev,h_RT_dev, (cudaStream_t) 0);
}


void mechanism_cuda::getIntEnergy_RT_mr_dev(const int nReactors, const double T_dev[], double u_RT_dev[]) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  static_cast<nasa_poly_group_cuda*>(thermo)->getH_RT_CUDA_mr(nReactors,T_dev,u_RT_dev, (cudaStream_t) 0);
  minusEqualOne(nSpc*nReactors,u_RT_dev);
}


void mechanism_cuda::getCfromVY_mr_dev(const int nReactors, const double v_dev[], const double y_dev[],
                               double c_dev[]) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  getCfromVY_mr_dev_wrapper(nReactors, nSpc, v_dev, y_dev, invMolWt_dev, c_dev);
}

void mechanism_cuda::getMassIntEnergyFromTY_mr_dev(const int nReactors, const double T_dev[], const double y_dev[],
                                                   double u_spc_dev[], double u_dev[] ) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  cudaMemsetAsync(u_dev,0,sizeof(double)*nReactors, (cudaStream_t) 0);
  static_cast<nasa_poly_group_cuda*>(thermo)->getH_RT_CUDA_mr(nReactors,T_dev, u_spc_dev, (cudaStream_t) 0);
  getMassIntEnergyFromTY_mr_dev_wrapper(nReactors, nSpc, T_dev, y_dev, RuInvMolWt_dev, u_spc_dev, u_dev);
}

void mechanism_cuda::getMassEnthalpyFromTY_mr_dev(const int nReactors, const double T_dev[], const double y_dev[],
                                                  double h_spc_dev[], double h_dev[] ) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  cudaMemsetAsync(h_dev,0,sizeof(double)*nReactors, (cudaStream_t) 0);
  static_cast<nasa_poly_group_cuda*>(thermo)->getH_RT_CUDA_mr(nReactors,T_dev, h_spc_dev, (cudaStream_t) 0);
  getMassEnthalpyFromTY_mr_dev_wrapper(nReactors, nSpc, T_dev, y_dev, RuInvMolWt_dev, h_spc_dev, h_dev);
}

void mechanism_cuda::getMolWtMixFromY_mr_dev(const int nReactors, const double y_dev[], double molWtMix_dev[]) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  cudaMemset(molWtMix_dev,0,sizeof(double)*nReactors);
  getMolWtMixFromY_mr_dev_wrapper(nReactors,nSpc,y_dev,invMolWt_dev,molWtMix_dev);
}

void mechanism_cuda::getDensityFromTPY_mr_dev(const int nReactors, const double T_dev[],
                                              const double P_dev[],const double y_dev[],
                                              double dens_dev[]) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  getMolWtMixFromY_mr_dev(nReactors,y_dev,dens_dev);
  getDensityfromTPMW_wrapper(nReactors,getGasConstant(),T_dev,P_dev,dens_dev);
}


void mechanism_cuda::getTemperatureFromEY_mr_dev(const int nReactors, const double E[], const double y[], double temperature_out[]) const
{
    assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
    const int maximum_iterations = 400;
    const double tolerance  = 1.0e-6;
    const double min_temperature = 100;
    const double max_temperature = 5000;

    thrust::device_vector<int> converged_flag(nReactors, 0);
    thrust::device_vector<double> temperature_min(nReactors, min_temperature);
    thrust::device_vector<double> temperature_max(nReactors, max_temperature);
    thrust::device_vector<double> tmp(nSpc*nReactors);
    thrust::device_vector<double> E1(nReactors);
    thrust::device_vector<double> E2(nReactors);
    thrust::device_vector<double> cv1(nReactors);
    thrust::device_vector<double> cv2(nReactors);
    getMassIntEnergyFromTY_mr_dev(nReactors, thrust::raw_pointer_cast(&temperature_min[0]), &y[0],
                                  thrust::raw_pointer_cast(&tmp[0]), thrust::raw_pointer_cast(&E1[0]));
    getMassCvFromTY_mr_dev(nReactors, thrust::raw_pointer_cast(&temperature_min[0]), &y[0],
                                      thrust::raw_pointer_cast(&tmp[0]),
                                      thrust::raw_pointer_cast(&cv1[0]));

    getMassIntEnergyFromTY_mr_dev(nReactors, thrust::raw_pointer_cast(&temperature_max[0]), &y[0],
                                  thrust::raw_pointer_cast(&tmp[0]), thrust::raw_pointer_cast(&E2[0]));
    getMassCvFromTY_mr_dev(nReactors, thrust::raw_pointer_cast(&temperature_max[0]), &y[0],
                                      thrust::raw_pointer_cast(&tmp[0]),
                                      thrust::raw_pointer_cast(&cv2[0]));

    getTemperatureFromEY_mr_dev_part1_wrapper(nReactors, min_temperature, max_temperature, &E[0],
                                              thrust::raw_pointer_cast(&E1[0]), 
                                              thrust::raw_pointer_cast(&E2[0]), 
                                              thrust::raw_pointer_cast(&cv1[0]), 
                                              thrust::raw_pointer_cast(&cv2[0]), 
                                              &temperature_out[0], 
                                              thrust::raw_pointer_cast(&converged_flag[0]));

    int sum = thrust::reduce(converged_flag.begin(), converged_flag.end(), (int) 0, thrust::plus<int>());
    if(sum == nReactors) return;

    for(int i = 0; i < maximum_iterations; ++i) {
        getMassIntEnergyFromTY_mr_dev(nReactors, &temperature_out[0], &y[0],
                                      thrust::raw_pointer_cast(&tmp[0]), thrust::raw_pointer_cast(&E1[0]));
        getMassCvFromTY_mr_dev(nReactors, &temperature_out[0], &y[0],
                               thrust::raw_pointer_cast(&tmp[0]),
                               thrust::raw_pointer_cast(&cv1[0]));
        getTemperatureFromEY_mr_dev_iter_wrapper(nReactors, tolerance, &E[0], thrust::raw_pointer_cast(&E1[0]),
                            thrust::raw_pointer_cast(&cv1[0]),
                            &temperature_out[0],
                            thrust::raw_pointer_cast(&converged_flag[0]));

        int sum = thrust::reduce(converged_flag.begin(), converged_flag.end(), (int) 0, thrust::plus<int>());
        if(sum == nReactors) break;
    }
}

void mechanism_cuda::getTemperatureFromHY_mr_dev(const int nReactors, const double H[], const double y[], double temperature_out[]) const
{
    assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
    const int maximum_iterations = 400;
    const double tolerance  = 1.0e-6;
    const double min_temperature = 100;
    const double max_temperature = 5000;

    thrust::device_vector<int> converged_flag(nReactors, 0);
    thrust::device_vector<double> temperature_min(nReactors, min_temperature);
    thrust::device_vector<double> temperature_max(nReactors, max_temperature);
    thrust::device_vector<double> tmp(nSpc*nReactors);
    thrust::device_vector<double> H1(nReactors);
    thrust::device_vector<double> H2(nReactors);
    thrust::device_vector<double> cp1(nReactors);
    thrust::device_vector<double> cp2(nReactors);
    getMassEnthalpyFromTY_mr_dev(nReactors, thrust::raw_pointer_cast(&temperature_min[0]), &y[0],
                                  thrust::raw_pointer_cast(&tmp[0]), thrust::raw_pointer_cast(&H1[0]));
    getMassCpFromTY_mr_dev(nReactors, thrust::raw_pointer_cast(&temperature_min[0]), &y[0],
                                      thrust::raw_pointer_cast(&tmp[0]),
                                      thrust::raw_pointer_cast(&cp1[0]));

    getMassEnthalpyFromTY_mr_dev(nReactors, thrust::raw_pointer_cast(&temperature_max[0]), &y[0],
                                  thrust::raw_pointer_cast(&tmp[0]), thrust::raw_pointer_cast(&H2[0]));
    getMassCpFromTY_mr_dev(nReactors, thrust::raw_pointer_cast(&temperature_max[0]), &y[0],
                                      thrust::raw_pointer_cast(&tmp[0]),
                                      thrust::raw_pointer_cast(&cp2[0]));

    getTemperatureFromEY_mr_dev_part1_wrapper(nReactors, min_temperature, max_temperature, &H[0],
                                              thrust::raw_pointer_cast(&H1[0]), 
                                              thrust::raw_pointer_cast(&H2[0]), 
                                              thrust::raw_pointer_cast(&cp1[0]), 
                                              thrust::raw_pointer_cast(&cp2[0]), 
                                              &temperature_out[0], 
                                              thrust::raw_pointer_cast(&converged_flag[0]));

    int sum = thrust::reduce(converged_flag.begin(), converged_flag.end(), (int) 0, thrust::plus<int>());
    if(sum == nReactors) return;

    for(int i = 0; i < maximum_iterations; ++i) {
        getMassEnthalpyFromTY_mr_dev(nReactors, &temperature_out[0], &y[0],
                                      thrust::raw_pointer_cast(&tmp[0]), thrust::raw_pointer_cast(&H1[0]));
        getMassCpFromTY_mr_dev(nReactors, &temperature_out[0], &y[0],
                               thrust::raw_pointer_cast(&tmp[0]),
                               thrust::raw_pointer_cast(&cp1[0]));
        getTemperatureFromEY_mr_dev_iter_wrapper(nReactors, tolerance, &H[0], thrust::raw_pointer_cast(&H1[0]),
                            thrust::raw_pointer_cast(&cp1[0]),
                            &temperature_out[0],
                            thrust::raw_pointer_cast(&converged_flag[0]));

        int sum = thrust::reduce(converged_flag.begin(), converged_flag.end(), (int) 0, thrust::plus<int>());
        if(sum == nReactors) break;
    }
}

void mechanism_cuda::getPressureFromTVY_mr_dev(const int nReactors, const double* T, const double* v,
                                               const double* y, double* P) const
{
  assert(nReactors <= static_cast<rate_const_cuda*>(Kconst)->nReactorsMax());
  getMolWtMixFromY_mr_dev(nReactors, y, P);
  getPressureFromTVY_wrapper(nReactors,getGasConstant(),T,v,P);
}


} // namespace zerork
