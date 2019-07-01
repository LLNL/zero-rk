#ifndef ZERORK_SPECIES_H
#define ZERORK_SPECIES_H

#include "element.h"

#include <vector>

namespace zerork {

class species
{
 public:
  species();
  species(const int idx, const char *nm, const std::vector<int> &count, 
          const std::vector<element> &list); 
  ~species();

  // data assignment
  void setSpecies(const int idx, const char *nm, const std::vector<int> &count,
                  const std::vector<element> &list); 

  // data access (constant member functions)
  int    getIndex() const {return index;}
  double getMolecularWeight() const {return molWt;}
  const char * getName_c_str() const {return name;}
  void printInfo() const;
  int getNumConstituents() const {return nConstituent;}
  int getTotalAtoms() const;
  int getCountOfAtomZ(const int zNum) const;

 private:
  int index;
  char *name;
  int nConstituent;
  double molWt;
  std::vector<int> constituentCount;
  std::vector<element> constituentList;
};

} // namespace zerork

#endif
