#include <stdio.h>
#include <stdlib.h> // needed for exit()
#include <string.h>
#include "element.h"
#include "atomicMassDB.h"
#include "utilities.h"

namespace zerork {

#ifdef EXIT_THROWS_EXCEPTION
  // create a local function to overide the system exit and throw an exception
  // with the status integer.
  static void exit(int status) {throw status;}
#endif // EXIT_THROWS_EXCEPTION

  
element::element()
{
  number=0;
  mass=0.0;
  strcpy(symbol,"**");
}

element::element(const int num, const double m, const char *sym)
{
  number = num;
  mass   = atomicMassDB[num-1];
  strncpy(symbol,sym,3);
  symbol[2]='\0'; // ensure that the symbol character array is terminated
}

// construct element with the weight and symbol defined in
// atomicMassDB.h
element::element(const int num)
{
  if(num >= 1 && num <= MAX_ATOMIC_NUM_DB)
    {
      number = num;
      mass   = atomicMassDB[num-1];
      strcpy(symbol,atomicSymbolDB[num-1]);
    }
  else if(num > MAX_ATOMIC_NUM_DB)
    {
      printf("ERROR: You should not be playing with element %d\n",num);
      exit(-1);   // replace hard quit with exception handling
    }
  else
    {
      printf("ERROR: element %d does not exist\n",num);
      exit(-1);   // replace hard quit with exception handling
    }
}
void element::setElement(const int num, const double m, const char *sym)
{
  number = num;
  mass   = atomicMassDB[num-1];
  strncpy(symbol,sym,3);
  symbol[2]='\0'; // ensure that the symbol character array is terminated
}

// set element using the atomic number with the weight and symbol defined
// in atomicMassDB.h
void element::setElement(const int num)
{
  if(num >= 1 && num <= MAX_ATOMIC_NUM_DB)
    {
      number = num;
      mass   = atomicMassDB[num-1];
      strcpy(symbol,atomicSymbolDB[num-1]);
    }
  else if(num > MAX_ATOMIC_NUM_DB)
    {
      printf("ERROR: You should not be playing with element %d\n",num);
      exit(-1);   // replace hard quit with exception handling
    }
  else
    {
      printf("ERROR: element %d does not exist\n",num);
      exit(-1);   // replace hard quit with exception handling
    }
}

// set element using the case insensitive symbol with the weight and atomic
// number defined in atomicMassDB.h
void element::setElement(const char *sym)
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
  if(foundNum==0)
    {
      printf("ERROR: could not find element %s (%s) in the database\n",
	     sym,lowerCaseSym);
      exit(-1);
    }
  
  number = foundNum;
  mass   = atomicMassDB[foundNum-1];
  strcpy(symbol,atomicSymbolDB[foundNum-1]);  
}

element::~element()
{}


void element::printInfo() const
{
  printf("# Element %2d  %2s  %10.6f kg/kmol\n",
	 number,symbol,mass);
}

} // namespace zerork

