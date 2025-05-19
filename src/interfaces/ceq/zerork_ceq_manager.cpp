
#include <memory>
#include <vector>

#include "zerork_ceq_manager.h"

namespace{
extern "C"
void eq_hp_init(int ne, int ns, double* mw, double* thermo, double* Ein);

extern "C"
void eq_hp(double p, double* T, double* y);

extern "C"
void eq_hp_free();
}

namespace zerork {

ceq_manager::ceq_manager(const zerork::mechanism& mech) {
    int num_species = mech.getNumSpecies();
    std::vector<double> thermoCoeffs(num_species*16);
    std::vector<double> thermoCoeffs_CEQ(num_species*15);
    mech.getThermoCoeffs(&thermoCoeffs[0]);
    for(int k=0; k<num_species; k++) {
      for(int l=0; l<15; l++) { //transpose
        thermoCoeffs_CEQ[l*num_species+k] = thermoCoeffs[k*16 + l];
      }
    }

    int ne = mech.getNumElements();
    std::vector<int> Ein_int(num_species*ne);
    std::vector<double> Ein_double(num_species*ne);
    mech.getSpeciesOxygenCount(&Ein_int[0]);
    mech.getSpeciesNitrogenCount(&Ein_int[2*num_species]);
    mech.getSpeciesHydrogenCount(&Ein_int[num_species]);
    if(ne>3) mech.getSpeciesCarbonCount(&Ein_int[3*num_species]);
    if(ne>4) mech.getSpeciesArgonCount(&Ein_int[4*num_species]);
    if(ne>5) mech.getSpeciesHeliumCount(&Ein_int[5*num_species]);

    for(int k=0; k<num_species*ne; k++) {
      Ein_double[k] = (double)Ein_int[k];
    }

    std::vector<double> molWt(num_species);
    mech.getMolWtSpc(&molWt[0]);

    //TODO: Thread safety in f90 lib (some handle to sys)
    eq_hp_init(ne, num_species, &molWt[0], &thermoCoeffs_CEQ[0], &Ein_double[0]);
}

ceq_manager::~ceq_manager() {
    eq_hp_free();
}

void ceq_manager::equilibrate(const double p, double& T, double* y) {
    eq_hp(p, &T, &y[0]);
}

} //namespace zerork
