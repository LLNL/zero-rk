
#include "multizone.h"

#include <algorithm>
#ifdef USE_MPI
#include "mpi.h"
#endif

namespace {
  int multizone_count = -1;
  std::map<int,int> zone_id_map;
  std::vector<int> reactor_zone_ids;
  std::vector<double> remap_factor;
} //anonymous namespace

void multizone_make_zones(const zerork::mechanism &mech,
                          int nreactors,
                          const double *T, const double *P,
                          const double *dpdt, const double *volume,
                          const double * massfracs, const double *cost,
                          const double * gpu,
                          int &nzones,
                          std::vector<double> &temp_zones,
                          std::vector<double> &press_zones,
                          std::vector<double> &dpdt_zones,
                          std::vector<double> &massfracs_zones,
                          std::vector<double> &cost_zones,
                          std::vector<double> &gpu_zones)
{
  //Future consideration.  Static zones
  multizone_count +=1;

  static int ncalls = 0;
  ncalls++;

  static int rank = 0;
  static int nranks = 1;
#ifdef USE_MPI
  if(ncalls == 1) {
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nranks);
  }
#endif

  int nSpc = mech.getNumSpecies();

  std::vector<double> t_profile;
  t_profile.push_back(0);
  t_profile.push_back(5000);
  std::vector<double> dt_prof;
  dt_prof.push_back(10);
  double phi_bin_size = 0.01;
  double T_cutoff = 0;

  int n_temp_zones=0;
  for(int i = 0; i < t_profile.size()-1; ++i)
  {
    n_temp_zones += (t_profile[i+1]-t_profile[i])/dt_prof[i];
  }
  n_temp_zones += 2; //Above and below t_profile

  reactor_zone_ids.resize(nreactors,-1);

  //Put each reactor in a zone
  int nzones_local = 0;
  for(int k = 0; k < nreactors; ++k) {
    int itemp = 0;
    int n_temp_ranges = t_profile.size()-1;
    double T_max_zone = t_profile[0];
    if( T[k] < t_profile[0] ) {itemp = 0;}
    else if (T[k] > t_profile[n_temp_ranges] ) { itemp = n_temp_zones; }
    else
    {
      int m = 0;
      itemp = 1;
      while( m < n_temp_ranges )
      {
        if( T[k] <= t_profile[m+1] )
        {
          int izone = (T[k]-t_profile[m])/dt_prof[m];
          T_max_zone = t_profile[m] + dt_prof[m]*(izone+1);
          itemp += izone;
          break;
        } else {
          itemp += (t_profile[m+1]-t_profile[m])/dt_prof[m];
        }
        m++;
      }
      if ( m == n_temp_ranges )
      {
          printf("ERROR: multizone_bins in temparature table for zones: \n"
                 "       Couldn't determine temperature zone for reactor %d "
                 " with temperature %f.\n",k,T[k]);
          exit(-1);
      }
    }
    if(T_max_zone <= T_cutoff)
    {
      //No zone for this cell
      reactor_zone_ids[k] = -1;
    }
    else
    {
      //Chemistry zones
      int ichem,ibig;
      bool phi_zoning = false;
      if(phi_zoning)
      {
        double phi_reactor = mech.getProgressEquivalenceRatio(&massfracs[k*nSpc]);
        printf("phi_reactor = %f\n",phi_reactor);
        ichem = phi_reactor/phi_bin_size;
        ibig = 0;
      }
      else
      {
        double smxfi,bmxfi;
        mech.getModifiedEquivalenceRatios(&massfracs[k*nSpc],smxfi,bmxfi);
        printf("smxf, bmxf= %f, %f\n",smxfi, bmxfi);
        ichem = smxfi/phi_bin_size;
        ibig = (bmxfi > 0.5) ? 0 : 1;
      }

      int zone_id = (n_temp_zones*ichem+itemp)*2+ibig;

      reactor_zone_ids[k] = zone_id;
    }
  }

  //Count the number of zones we have locally
  std::vector<int> unique_zone_ids_local = reactor_zone_ids;
  std::sort(unique_zone_ids_local.begin(), unique_zone_ids_local.end());
  std::vector<int>::iterator last_unique =
      std::unique(unique_zone_ids_local.begin(),unique_zone_ids_local.end());
  unique_zone_ids_local.erase(last_unique, unique_zone_ids_local.end());

  nzones_local = unique_zone_ids_local.size();

#ifdef USE_MPI
  //Communicate zones and generate global list.
  std::vector<int> nzones_ranks;
  if(rank == 0) nzones_ranks.resize(nranks);
  MPI_Gather(&nzones_local,1,MPI_INT,&nzones_ranks[0],1,MPI_INT,0,MPI_COMM_WORLD);

  std::vector<int> unique_zone_ids;
  std::vector<int> recv_displs;
  nzones = 0;
  if(rank == 0) {
    recv_displs.resize(nranks);
    for(int i = 0; i < nranks; ++i) {
      nzones += nzones_ranks[i];
      recv_displs[i] = nzones - nzones_ranks[i];
    }
    unique_zone_ids.resize(nzones);
  }
  MPI_Gatherv(&unique_zone_ids_local[0], nzones_local, MPI_INT,
              &unique_zone_ids[0], &nzones_ranks[0],
              &recv_displs[0],MPI_INT,0,MPI_COMM_WORLD);
  if(rank == 0) {
    std::sort(unique_zone_ids.begin(),unique_zone_ids.end());
    last_unique = std::unique(unique_zone_ids.begin(), unique_zone_ids.end());
    unique_zone_ids.erase(last_unique,unique_zone_ids.end());
    nzones = unique_zone_ids.size();
  }
  MPI_Bcast(&nzones,1,MPI_INT,0,MPI_COMM_WORLD);
  unique_zone_ids.resize(nzones);
  MPI_Bcast(&unique_zone_ids[0],nzones,MPI_INT,0,MPI_COMM_WORLD);
#else
  std::vector<int> &unique_zone_ids = unique_zone_ids_local;
  nzones = nzones_local;
#endif

  zone_id_map.clear();

  for(int i = 0; i < nzones; ++i)
  {
    zone_id_map[unique_zone_ids[i]] = i;
  }

  zone_id_map[-1] = -1;


  //Passed by reference.
  temp_zones.resize(nzones,0.0);
  press_zones.resize(nzones,0.0);
  dpdt_zones.resize(nzones,0.0);
  massfracs_zones.resize(nzones*nSpc,0.0);
  cost_zones.resize(nzones,0.0);
  gpu_zones.resize(nzones,0.0);


  std::vector<double> mass_zones_local(nzones,0.0);
  std::vector<double> Cv_zones_local(nzones,0.0);
  std::vector<double> volume_zones_local(nzones,0.0);
  std::vector<int> count_zones_local(nzones,0);

#ifdef USE_MPI
  std::vector<double> temp_zones_local(nzones,0.0);
  std::vector<double> press_zones_local(nzones,0.0);
  std::vector<double> dpdt_zones_local(nzones,0.0);
  std::vector<double> massfracs_zones_local(nzones*nSpc,0.0);
  std::vector<double> cost_zones_local(nzones,0.0);
  std::vector<double> gpu_zones_local(nzones,0.0);

  std::vector<double> mass_zones(nzones,0.0);
  std::vector<double> Cv_zones(nzones,0.0);
  std::vector<double> volume_zones(nzones,0.0);
  std::vector<int> count_zones(nzones,0);
#else
  std::vector<double> &temp_zones_local = temp_zones;
  std::vector<double> &press_zones_local = press_zones;
  std::vector<double> &dpdt_zones_local = dpdt_zones;
  std::vector<double> &massfracs_zones_local = massfracs_zones;
  std::vector<double> &cost_zones_local = cost_zones;
  std::vector<double> &gpu_zones_local = gpu_zones;

  std::vector<double> &mass_zones = mass_zones_local;
  std::vector<double> &volume_zones = volume_zones_local;
  std::vector<double> &Cv_zones = Cv_zones_local;
  std::vector<int> &count_zones = count_zones_local;
#endif

  for(int k = 0; k < nreactors; ++k)
  {
    int izone = zone_id_map[reactor_zone_ids[k]];
    if(izone != -1) {
      double reactor_density = mech.getDensityFromTPY(T[k],P[k],
                                                      &massfracs[k*nSpc]);
      double reactor_Cv = mech.getMassCvFromTY(T[k],&massfracs[k*nSpc]);
      double reactor_volume = volume[k];
      double reactor_mass = reactor_density*reactor_volume;

      volume_zones_local[izone] += reactor_volume;
      mass_zones_local[izone] += reactor_mass;
      Cv_zones_local[izone] += reactor_Cv*reactor_mass;
      count_zones_local[izone] += 1;

      temp_zones_local[izone] += T[k]*reactor_Cv*reactor_mass;
      press_zones_local[izone] += P[k]*reactor_volume;
      dpdt_zones_local[izone] += dpdt[k]*reactor_volume;
      for(int i = 0; i < nSpc; ++i)
      {
        massfracs_zones_local[izone*nSpc+i] += massfracs[k*nSpc+i]*reactor_mass;
      }
      cost_zones_local[izone] += cost[k];
      gpu_zones_local[izone] += gpu[k];
    }
  }

#ifdef USE_MPI
  //Communicate zone info
  MPI_Allreduce(&volume_zones_local[0],&volume_zones[0],nzones,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&mass_zones_local[0],&mass_zones[0],nzones,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&Cv_zones_local[0],&Cv_zones[0],nzones,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&temp_zones_local[0],&temp_zones[0],nzones,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&press_zones_local[0],&press_zones[0],nzones,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&dpdt_zones_local[0],&dpdt_zones[0],nzones,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&massfracs_zones_local[0],&massfracs_zones[0],nzones*nSpc,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&cost_zones_local[0],&cost_zones[0],nzones,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&gpu_zones_local[0],&gpu_zones[0],nzones,
                MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&count_zones_local[0],&count_zones[0],nzones,
                MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif

  //Normalize after communication
  for(int izone = 0; izone < nzones; ++izone)
  {
    temp_zones[izone] /= Cv_zones[izone]; //Cv_zones includes zone mass
    press_zones[izone] /= volume_zones[izone];
    dpdt_zones[izone] /= volume_zones[izone];
    for(int i = 0; i < nSpc; ++i) {
      massfracs_zones[izone*nSpc+i] /= mass_zones[izone];
    }
    cost_zones[izone] /= count_zones[izone];
    gpu_zones[izone] /= count_zones[izone];
  }

  //Calc remap factor
  remap_factor.resize(nreactors,0.0);
  for(int k = 0; k < nreactors; ++k)
  {
    int izone = zone_id_map[reactor_zone_ids[k]];
    if(izone != -1) {
      double ch_reactor = mech.getCHValue(&massfracs[nSpc*k]);
      double ch_zone = mech.getCHValue(&massfracs_zones[nSpc*izone]);
      remap_factor[k] = ch_reactor / ch_zone;
    } else {
      remap_factor[k] = 1.0;
    }
  }
}


void multizone_remap(const zerork::mechanism &mech,
                     int nzones,
                     const std::vector<double> &massfracs_zones0,
                     const std::vector<double> &massfracs_zones,
                     int nreactors,
                     double* massfracs,
                     const std::vector<double> &cost_zones,
                     double* cost,
                     const std::vector<double> &gpu_zones,
                     double* gpu)
{
  static int num_remaps = -1;
  num_remaps += 1;
  if(multizone_count != num_remaps)
  {
    printf("ERROR: multizone_remap called inconsistently.\n");
    exit(-1);
  }

  int nSpc = mech.getNumSpecies();
  //Simplistic remap to start.
  for(int k = 0; k < nreactors; ++k)
  {
    int izone = zone_id_map[reactor_zone_ids[k]];
    if(izone != -1)
    {
      double ytotal = 0.0;
      for (int i=0; i<nSpc; i++)
      {
        double dyi=massfracs_zones[nSpc*izone+i]-massfracs_zones0[nSpc*izone+i];
        if(dyi < 0.0)
        {
          if(massfracs_zones0[nSpc*izone+i] > 0.0)
          {
            double dyi_reactor = dyi*massfracs[nSpc*k+i]/massfracs_zones0[nSpc*izone+i];
            massfracs[nSpc*k+i] += dyi_reactor;
            if(massfracs[nSpc*k+i] < 0.0) massfracs[nSpc*k+i] = 0.0;
          }
        }
        else
        {
          massfracs[nSpc*k+i] += dyi;
        }
        ytotal += massfracs[nSpc*k+i];
      }
      for (int i=0; i<nSpc; i++)
      {
        massfracs[nSpc*k+i] /= ytotal;
      }
      //Update cost & gpu
      cost[k] = cost_zones[izone];
      gpu[k] = gpu_zones[izone];
    }
  }
}

void multizone_remap_full(const zerork::mechanism &mech,
                     int nzones,
                     const std::vector<double> &massfracs_zones0,
                     const std::vector<double> &massfracs_zones,
                     int nreactors,
                     double* massfracs,
                     const std::vector<double> &cost_zones,
                     double* cost,
                     const std::vector<double> &gpu_zones,
                     double* gpu)
{
  static int num_remaps = -1;
  num_remaps += 1;
  if(multizone_count != num_remaps)
  {
    printf("ERROR: multizone_remap called inconsistently.\n");
    exit(-1);
  }

  int nSpc = mech.getNumSpecies();

  int idx_o2 = -1;
  int idx_co2 = -1;
  int idx_h2o = -1;
  int idx_n2 = -1;
  std::vector<int> numC(nSpc);
  std::vector<int> numH(nSpc);
  std::vector<int> numO(nSpc);
  std::vector<double> molWt(nSpc);

  mech.getCarbonAtomCount(&numC[0]);
  mech.getHydrogenAtomCount(&numH[0]);
  mech.getOxygenAtomCount(&numO[0]);
  mech.getMolWtSpc(&molWt[0]);

  for(int i = 0; i < nSpc; ++i)
  {
      if(strcmp("O2",mech.getSpeciesName(i))==0 ||
         strcmp("o2",mech.getSpeciesName(i))==0 )
         {
           idx_o2 = i;
         }
      else if(strcmp("CO2",mech.getSpeciesName(i))==0 ||
              strcmp("co2",mech.getSpeciesName(i))==0 )
         {
           idx_co2 = i;
         }
      else if(strcmp("H2O",mech.getSpeciesName(i))==0 ||
              strcmp("h2o",mech.getSpeciesName(i))==0 )
         {
           idx_h2o = i;
         }
      else if(strcmp("N2",mech.getSpeciesName(i))==0 ||
              strcmp("n2",mech.getSpeciesName(i))==0 )
         {
           idx_n2 = i;
         }
  }
  if(idx_o2 < 0 || idx_co2 < 0 || idx_h2o < 0 || idx_n2 < 0)
  {
    printf("ERROR: Failed to find key species %d %d %d %d\n",idx_o2,idx_co2,idx_h2o,idx_n2);
    exit(-1);
  }

  std::vector<int> zone_flags_local(nzones,0);
  std::vector<double> massfracs0(nSpc);
  for(int k = 0; k < nreactors; ++k)
  {
    int izone = zone_id_map[reactor_zone_ids[k]];
    if(izone != -1 && zone_flags_local[izone] == 0)
    {
      double tol = 1.0e-10;
      memcpy(&massfracs0[0],&massfracs[k*nSpc],nSpc*sizeof(double));
      double ytotal = 0.0;
      double totC_old = 0.0;
      double totH_old = 0.0;
      double totO_old = 0.0;
      double totC_new = 0.0;
      double totH_new = 0.0;
      double totO_new = 0.0;
      double fact = remap_factor[k];
      for (int i=0; i<nSpc; i++)
      {
        double invMolWt = 1.0/molWt[i];
        totC_old += massfracs[nSpc*k+i]*numC[i]*invMolWt;
        totH_old += massfracs[nSpc*k+i]*numH[i]*invMolWt;
        totO_old += massfracs[nSpc*k+i]*numO[i]*invMolWt;
        if(i != idx_o2 && i != idx_co2 && i != idx_h2o && i != idx_n2)
        {
          //massfracs[nSpc*k+i] = max(fact * massfracs_zones[nSpc*izone+i],0.0);
          massfracs[nSpc*k+i] = fact*massfracs_zones[nSpc*izone+i];
          totC_new += massfracs[nSpc*k+i]*numC[i]*invMolWt;
          totH_new += massfracs[nSpc*k+i]*numH[i]*invMolWt;
          totO_new += massfracs[nSpc*k+i]*numO[i]*invMolWt;
          ytotal += massfracs[nSpc*k+i];
        }
      }
      massfracs[nSpc*k+idx_co2] = (totC_old-totC_new)/numC[idx_co2]*molWt[idx_co2];
      if(massfracs[nSpc*k+idx_co2] < -tol) {
        zone_flags_local[izone] = 1;
      }
      massfracs[nSpc*k+idx_co2] = max(0.0,massfracs[nSpc*k+idx_co2]);

      massfracs[nSpc*k+idx_h2o] = (totH_old-totH_new)/numH[idx_h2o]*molWt[idx_h2o];
      if(massfracs[nSpc*k+idx_h2o] < -tol) {
        zone_flags_local[izone] = 1;
      }
      massfracs[nSpc*k+idx_h2o] = max(0.0,massfracs[nSpc*k+idx_h2o]);

      totO_new += massfracs[nSpc*k+idx_co2]*(numO[idx_co2]/molWt[idx_co2]);
      totO_new += massfracs[nSpc*k+idx_h2o]*(numO[idx_h2o]/molWt[idx_h2o]);
      massfracs[nSpc*k+idx_o2] = (totO_old-totO_new)/numO[idx_o2]*molWt[idx_o2];
      if(massfracs[nSpc*k+idx_o2] < -tol) {
        zone_flags_local[izone] = 1;
      }
      massfracs[nSpc*k+idx_o2] = max(0.0,massfracs[nSpc*k+idx_o2]);

      ytotal += massfracs[nSpc*k+idx_co2];
      ytotal += massfracs[nSpc*k+idx_h2o];
      ytotal += massfracs[nSpc*k+idx_o2];
      massfracs[nSpc*k+idx_n2] = 1.0 - ytotal;
      if(massfracs[nSpc*k+idx_n2] < -tol) {
        zone_flags_local[izone] = 1;
      }
      massfracs[nSpc*k+idx_n2] = max(0.0,massfracs[nSpc*k+idx_n2]);

      //Set back to initial mass fracs and re-do remap in backup mode
      if(zone_flags_local[izone] ==1) {
        memcpy(&massfracs[k*nSpc],&massfracs0[0],nSpc*sizeof(double));
      }
      //Update cost & gpu
      cost[k] = cost_zones[izone];
      gpu[k] = gpu_zones[izone];
    }
  }

#ifdef USE_MPI
  std::vector<int> zone_flags(nzones,0);
  MPI_Allreduce(&zone_flags_local[0],&zone_flags[0],nzones,
                MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#else
  std::vector<int> &zone_flags = zone_flags_local;
#endif

#if 1
  for(int k = 0; k < nreactors; ++k)
  {
    int izone = zone_id_map[reactor_zone_ids[k]];
    if(izone != -1 && zone_flags[izone] == 1)
    {
      printf("MULTIZONE RECOVERING\n");
      double ytotal = 0.0;
      for (int i=0; i<nSpc; i++)
      {
        double dyi=massfracs_zones[nSpc*izone+i]-massfracs_zones0[nSpc*izone+i];
        if(dyi < 0.0)
        {
          if(massfracs_zones0[nSpc*izone+i] > 0.0)
          {
            double dyi_reactor = dyi*massfracs[k*nSpc+i]/massfracs_zones0[nSpc*izone+i];
            massfracs[nSpc*k+i] += dyi_reactor;
            if(massfracs[nSpc*k+i] < 0.0) massfracs[nSpc*k+i] = 0.0;
          }
        }
        else
        {
          massfracs[nSpc*k+i] += dyi;
        }
        ytotal += massfracs[nSpc*k+i];
      }
      for (int i=0; i<nSpc; i++)
      {
        massfracs[nSpc*k+i] /= ytotal;
      }
    }
  }
#endif
}

void multizone_remap_full2(const zerork::mechanism &mech,
                     int nzones,
                     const std::vector<double> &massfracs_zones0,
                     const std::vector<double> &massfracs_zones,
                     int nreactors,
                     double* massfracs,
                     const std::vector<double> &cost_zones,
                     double* cost,
                     const std::vector<double> &gpu_zones,
                     double* gpu)
{
  static int num_remaps = -1;
  num_remaps += 1;
  if(multizone_count != num_remaps)
  {
    printf("ERROR: multizone_remap called inconsistently.\n");
    exit(-1);
  }

  int nSpc = mech.getNumSpecies();

  int idx_o2 = -1;
  int idx_co2 = -1;
  int idx_h2o = -1;
  int idx_n2 = -1;
  std::vector<int> numC(nSpc);
  std::vector<int> numH(nSpc);
  std::vector<int> numO(nSpc);
  std::vector<double> molWt(nSpc);

  mech.getCarbonAtomCount(&numC[0]);
  mech.getHydrogenAtomCount(&numH[0]);
  mech.getOxygenAtomCount(&numO[0]);
  mech.getMolWtSpc(&molWt[0]);

  for(int i = 0; i < nSpc; ++i)
  {
      if(strcmp("O2",mech.getSpeciesName(i))==0 ||
         strcmp("o2",mech.getSpeciesName(i))==0 )
         {
           idx_o2 = i;
         }
      else if(strcmp("CO2",mech.getSpeciesName(i))==0 ||
              strcmp("co2",mech.getSpeciesName(i))==0 )
         {
           idx_co2 = i;
         }
      else if(strcmp("H2O",mech.getSpeciesName(i))==0 ||
              strcmp("h2o",mech.getSpeciesName(i))==0 )
         {
           idx_h2o = i;
         }
      else if(strcmp("N2",mech.getSpeciesName(i))==0 ||
              strcmp("n2",mech.getSpeciesName(i))==0 )
         {
           idx_n2 = i;
         }
  }
  if(idx_o2 < 0 || idx_co2 < 0 || idx_h2o < 0 || idx_n2 < 0)
  {
    printf("ERROR: Failed to find key species %d %d %d %d\n",idx_o2,idx_co2,idx_h2o,idx_n2);
    exit(-1);
  }

  std::vector<int> zone_flags_local(nzones,0);
  std::vector<double> massfracs0(nSpc);
  for(int k = 0; k < nreactors; ++k)
  {
    int izone = zone_id_map[reactor_zone_ids[k]];
    if(izone != -1 && zone_flags_local[izone] == 0)
    {
      double tol = 1.0e-10;
      memcpy(&massfracs0[0],&massfracs[k*nSpc],nSpc*sizeof(double));
      double ytotal = 0.0;
      double totC_old = 0.0;
      double totH_old = 0.0;
      double totO_old = 0.0;
      double totC_new = 0.0;
      double totH_new = 0.0;
      double totO_new = 0.0;
      double fact = remap_factor[k];
      for (int i=0; i<nSpc; i++)
      {
        double invMolWt = 1.0/molWt[i];
        totC_old += massfracs[nSpc*k+i]*numC[i]*invMolWt;
        totH_old += massfracs[nSpc*k+i]*numH[i]*invMolWt;
        totO_old += massfracs[nSpc*k+i]*numO[i]*invMolWt;
        if(i != idx_o2 && i != idx_co2 && i != idx_h2o && i != idx_n2)
        {
          double dyi=massfracs_zones[nSpc*izone+i]-massfracs_zones0[nSpc*izone+i];
          if(dyi < 0.0)
          {
            if(massfracs_zones0[nSpc*izone+i] > 0.0)
            {
              double dyi_reactor = dyi*massfracs[k*nSpc+i]/(massfracs_zones0[nSpc*izone+i]+1.0e-30);
              massfracs[nSpc*k+i] += dyi_reactor;
              if(massfracs[nSpc*k+i] < 0.0) massfracs[nSpc*k+i] = 0.0;
            }
          }
          else
          {
            massfracs[nSpc*k+i] += dyi*fact;
          }
          totC_new += massfracs[nSpc*k+i]*numC[i]*invMolWt;
          totH_new += massfracs[nSpc*k+i]*numH[i]*invMolWt;
          totO_new += massfracs[nSpc*k+i]*numO[i]*invMolWt;
          ytotal += massfracs[nSpc*k+i];
        }
      }
      massfracs[nSpc*k+idx_co2] = (totC_old-totC_new)/numC[idx_co2]*molWt[idx_co2];
      if(massfracs[nSpc*k+idx_co2] < -tol) {
        zone_flags_local[izone] = 1;
      }
      massfracs[nSpc*k+idx_co2] = max(0.0,massfracs[nSpc*k+idx_co2]);

      massfracs[nSpc*k+idx_h2o] = (totH_old-totH_new)/numH[idx_h2o]*molWt[idx_h2o];
      if(massfracs[nSpc*k+idx_h2o] < -tol) {
        zone_flags_local[izone] = 1;
      }
      massfracs[nSpc*k+idx_h2o] = max(0.0,massfracs[nSpc*k+idx_h2o]);

      totO_new += massfracs[nSpc*k+idx_co2]*(numO[idx_co2]/molWt[idx_co2]);
      totO_new += massfracs[nSpc*k+idx_h2o]*(numO[idx_h2o]/molWt[idx_h2o]);
      massfracs[nSpc*k+idx_o2] = (totO_old-totO_new)/numO[idx_o2]*molWt[idx_o2];
      if(massfracs[nSpc*k+idx_o2] < -tol) {
        zone_flags_local[izone] = 1;
      }
      massfracs[nSpc*k+idx_o2] = max(0.0,massfracs[nSpc*k+idx_o2]);

      ytotal += massfracs[nSpc*k+idx_co2];
      ytotal += massfracs[nSpc*k+idx_h2o];
      ytotal += massfracs[nSpc*k+idx_o2];
      massfracs[nSpc*k+idx_n2] = 1.0 - ytotal;
      if(massfracs[nSpc*k+idx_n2] < -tol) {
        zone_flags_local[izone] = 1;
      }
      massfracs[nSpc*k+idx_n2] = max(0.0,massfracs[nSpc*k+idx_n2]);
      //Set back to initial mass fracs and re-do remap in backup mode
      if(zone_flags_local[izone] == 1) {
        memcpy(&massfracs[k*nSpc],&massfracs0[0],nSpc*sizeof(double));
      }
      //Update cost & gpu
      cost[k] = cost_zones[izone];
      gpu[k] = gpu_zones[izone];
    }
  }

#ifdef USE_MPI
  std::vector<int> zone_flags(nzones,0);
  MPI_Allreduce(&zone_flags_local[0],&zone_flags[0],nzones,
                MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#else
  std::vector<int> &zone_flags = zone_flags_local;
#endif

#if 1
  for(int k = 0; k < nreactors; ++k)
  {
    int izone = zone_id_map[reactor_zone_ids[k]];
    if(izone != -1 && zone_flags[izone] == 1)
    {
      printf("MULTIZONE RECOVERING\n");
      double ytotal = 0.0;
      for (int i=0; i<nSpc; i++)
      {
        double dyi=massfracs_zones[nSpc*izone+i]-massfracs_zones0[nSpc*izone+i];
        if(dyi < 0.0)
        {
          if(massfracs_zones0[nSpc*izone+i] > 0.0)
          {
            double dyi_reactor = dyi*massfracs[k*nSpc+i]/massfracs_zones0[nSpc*izone+i];
            massfracs[nSpc*k+i] += dyi_reactor;
            if(massfracs[nSpc*k+i] < 0.0) massfracs[nSpc*k+i] = 0.0;
          }
        }
        else
        {
          massfracs[nSpc*k+i] += dyi;
        }
        ytotal += massfracs[nSpc*k+i];
      }
      for (int i=0; i<nSpc; i++)
      {
        massfracs[nSpc*k+i] /= ytotal;
      }
    }
  }
#endif
}
