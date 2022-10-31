#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "element.h"
#include "atomicMassDB.h"
#include "utilities.h"

namespace zerork {
  
element::element()
{
  number=0;
  mass=0.0;
  strcpy(symbol,"**");
}

element::element(const int num, const double m, const char *sym)
{
  assert(("Element out of range", num >=1 && num <= MAX_ATOMIC_NUM_DB));
  number = num;
  mass   = atomicMassDB[num-1];
  strncpy(symbol,sym,3);
  symbol[2]='\0'; // ensure that the symbol character array is terminated
}

// construct element with the weight and symbol defined in
// atomicMassDB.h
element::element(const int num)
{
  assert(("Element out of range", num >=1 && num <= MAX_ATOMIC_NUM_DB));
  number = num;
  mass   = atomicMassDB[num-1];
  strcpy(symbol,atomicSymbolDB[num-1]);
}

bool element::setElement(const int num, const double m, const char *sym)
{
  assert(("Element out of range", num >=1 && num <= MAX_ATOMIC_NUM_DB));
  number = num;
  mass   = atomicMassDB[num-1];
  strncpy(symbol,sym,3);
  symbol[2]='\0'; // ensure that the symbol character array is terminated
  return true;
}

// set element using the atomic number with the weight and symbol defined
// in atomicMassDB.h
bool element::setElement(const int num)
{
  if(num < 1 && num > MAX_ATOMIC_NUM_DB) {
    return false;
  }
  number = num;
  mass   = atomicMassDB[num-1];
  strcpy(symbol,atomicSymbolDB[num-1]);
  return true;
}

// set element using the case insensitive symbol with the weight and atomic
// number defined in atomicMassDB.h
bool element::setElement(const char *sym)
{
  int j,foundNum;
  char lowerCaseSym[3],lowerCaseDB[3];
  
  strncpy(lowerCaseSym,sym,3);
  lowerCaseSym[2]='\0';
  upperToLower(lowerCaseSym);
  
  foundNum=0;
  j=0;
  while(!foundNum && j<MAX_ATOMIC_NUM_DB)
    {
      strcpy(lowerCaseDB,atomicSymbolDB[j]);
      upperToLower(lowerCaseDB);
      if(strncmp(lowerCaseDB,lowerCaseSym,2)==0)
	{foundNum=j+1;}
      j++;
    }
  if(foundNum==0) {
    return false;
  }
  
  number = foundNum;
  mass   = atomicMassDB[foundNum-1];
  strcpy(symbol,atomicSymbolDB[foundNum-1]);  
  return true;
}

element::~element()
{}


void element::printInfo() const
{
  printf("# Element %2d  %2s  %10.6f kg/kmol\n",
	 number,symbol,mass);
}

} // namespace zerork

