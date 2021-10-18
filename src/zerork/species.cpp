#include "species.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

namespace zerork {

species::species()
{
  index=-1;
  nConstituent=1;
  molWt=0.0;
  constituentCount.resize(1);
  constituentCount[0]=1;
  constituentList.resize(1);
  constituentList[0]=element();
  name=(char *)malloc(sizeof(char)*10);
  strcpy(name,"undefined");
}

species::species(const int idx, const char *nm,
                 const std::vector<int> &count, 
                 const std::vector<element> &list)
{
  int j;
  int nameLength=strlen(nm)+1;

  if(name!=NULL)
    {free(name);}
  name=(char *)malloc(sizeof(char)*nameLength);
  strcpy(name,nm);
  index=idx;
  nConstituent=count.size();
  constituentCount.resize(nConstituent);
  constituentList.resize(nConstituent);

  molWt=0.0;

  for(j=0; j<nConstituent; j++)
  {
      constituentCount[j]=count[j];
      constituentList[j]=list[j];
      molWt+=constituentCount[j]*constituentList[j].getMass();
  }

}


species::~species()
{free(name);}

void species::setSpecies(const int idx, const char *nm, 
                         const std::vector<int> &count, 
                         const std::vector<element> &list)
{
  int j;
  int nameLength=strlen(nm)+1;

  if(name!=NULL)
    {free(name);}
  name=(char *)malloc(sizeof(char)*nameLength);
  strncpy(name,nm,nameLength);
  index=idx;
  nConstituent=count.size();
  constituentCount.resize(nConstituent);
  constituentList.resize(nConstituent);

  molWt=0.0;

  for(j=0; j<nConstituent; j++)
  {
      constituentCount[j]=count[j];
      constituentList[j]=list[j];
      molWt+=constituentCount[j]*constituentList[j].getMass();
  }

}

void species::printInfo() const
{
  int j;
  printf("# Species %d %s %10.6f kg/kmol\n",index,name,molWt);
  for(j=0; j<nConstituent; j++)
  {
      printf("#   composed of %d*%s\n",constituentCount[j],
             constituentList[j].getSymbol_c_str());
  }
}

int species::getCountOfAtomZ(const int zNum) const
{
  int j;
  int atomCtr=0;
  for(j=0; j<nConstituent; j++)
  {
      if(constituentList[j].getNumber() == zNum)
        {atomCtr+=constituentCount[j];}
  }
  return atomCtr;
}

int  species::getTotalAtoms() const
{
  int j;
  int atomCtr=0;
  for(j=0; j<nConstituent; j++)
    {atomCtr+=constituentCount[j];}

  return atomCtr;
}

} // namespace zerork

