
#include <memory>
#include <vector>

#include "zerork_cantera_manager.h"

#include "cantera/thermo/Species.h"
#include "cantera/thermo/NasaPoly2.h"

static int ns_cant;

namespace zerork {

cantera_manager::cantera_manager(const zerork::mechanism& mech) {

  int num_species = mech.getNumSpecies();

  std::vector<double> thermoCoef(num_species*LDA_THERMO_POLY_D5R2);
  mech.getThermoCoeffs(&thermoCoef[0]);
  std::map<std::string, std::map<std::string, int>> speciesElemInfo = mech.getSpeciesElementInfo();
  for(int k=0; k<num_species; k++) {
    const char * name = mech.getSpeciesName(k);
    std::map<std::string, int> elemMap = speciesElemInfo[name];
    std::map<std::string, double> canteraElemMap;
    for (auto const& e : elemMap) {
      canteraElemMap[e.first] = (double) e.second;
    }
    
    auto s = make_shared<Cantera::Species>(name, canteraElemMap);
    //const double* coefs = &(thermoCoef[k*LDA_THERMO_POLY_D5R2]);
    std::vector<double> coefs(15);
    coefs[0] = thermoCoef[k*LDA_THERMO_POLY_D5R2+0];
    coefs[1] = thermoCoef[k*LDA_THERMO_POLY_D5R2+8];
    coefs[2] = thermoCoef[k*LDA_THERMO_POLY_D5R2+9];
    coefs[3] = thermoCoef[k*LDA_THERMO_POLY_D5R2+10];
    coefs[4] = thermoCoef[k*LDA_THERMO_POLY_D5R2+11];
    coefs[5] = thermoCoef[k*LDA_THERMO_POLY_D5R2+12];
    coefs[6] = thermoCoef[k*LDA_THERMO_POLY_D5R2+13];
    coefs[7] = thermoCoef[k*LDA_THERMO_POLY_D5R2+14];
    coefs[8] = thermoCoef[k*LDA_THERMO_POLY_D5R2+1];
    coefs[9] = thermoCoef[k*LDA_THERMO_POLY_D5R2+2];
    coefs[10] = thermoCoef[k*LDA_THERMO_POLY_D5R2+3];
    coefs[11] = thermoCoef[k*LDA_THERMO_POLY_D5R2+4];
    coefs[12] = thermoCoef[k*LDA_THERMO_POLY_D5R2+5];
    coefs[13] = thermoCoef[k*LDA_THERMO_POLY_D5R2+6];
    coefs[14] = thermoCoef[k*LDA_THERMO_POLY_D5R2+7];
    
    s->thermo = make_shared<Cantera::NasaPoly2>(200, 3500, 101325, &coefs[0]);
    gas.addSpecies(s);
  }
  gas.initThermo();
  ns_cant = num_species;
}

cantera_manager::~cantera_manager() {}

void cantera_manager::equilibrate(const double p, double& T, double* y) {
  gas.setState_TPY(T, p, y);
  gas.getMassFractions(y);
  printf("p: %g\n", p);
  printf("Tpre: %g\n", T);
  for(int i = 0; i < ns_cant; ++i) {
    printf("%g\t", y[i]);
  }
  printf("\n");
  gas.equilibrate("HP", "element_potential");
  //gas.equilibrate("HP", "vcs");
  gas.getMassFractions(y);
  T = gas.temperature();
  printf("Teq: %g\n", T);
  for(int i = 0; i < ns_cant; ++i) {
    printf("%g\t", y[i]);
  }
  printf("\n");
}

} //namespace zerork
