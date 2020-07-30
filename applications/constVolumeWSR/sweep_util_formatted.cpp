#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sweep_util_formatted.h"

idt_sweep_params::idt_sweep_params(char *inputFileName)
{
  int j,k;
  FILE *inputFptr=fopen(inputFileName,"r");

  if(inputFptr == NULL)
    {
      printf("ERROR: could not open input file %s for read\n",inputFileName);
      fflush(stdout); exit(-1);
    }

  // read in the filenames
  mechFileName  = readFileLine_first_string(inputFptr);
  thermFileName = readFileLine_first_string(inputFptr);
  logFileName   = readFileLine_first_string(inputFptr);
  idtFileName   = readFileLine_first_string(inputFptr);
  thistFileName = readFileLine_first_string(inputFptr);

  gasMech = new zerork::mechanism(mechFileName,thermFileName,logFileName);

  // ----- parse the fuel and oxygen mixtures -----
  nSpc=gasMech->getNumSpecies();
  fuelMoleFrac = new double[nSpc];
  oxidMoleFrac = new double[nSpc];
  traceMoleFrac = new double[nSpc];

  // set the initial mole fractions to a large negative number to identify
  // species declared with zero mole fraction for tracking purposes
  for(j=0; j<nSpc; j++)
    {fuelMoleFrac[j]=oxidMoleFrac[j]=traceMoleFrac[j]=-1.0e300;} // use this to identify
  
  nFuelSpc=readFileLine_first_int(inputFptr);
  fuelSpcIdx = new int[nFuelSpc];
  for(j=0; j<nFuelSpc; j++)
    {readFileLine_first_comp(inputFptr,*gasMech,fuelMoleFrac);}
  nOxidSpc=readFileLine_first_int(inputFptr);
  oxidSpcIdx = new int[nOxidSpc];
  for(j=0; j<nOxidSpc; j++)
    {readFileLine_first_comp(inputFptr,*gasMech,oxidMoleFrac);}
  nTraceSpc=readFileLine_first_int(inputFptr);
  traceSpcIdx = new int[nTraceSpc];
  for(j=0; j<nTraceSpc; j++)
    {readFileLine_first_comp(inputFptr,*gasMech,traceMoleFrac);}

  // record the species indexes found for the fuel
  k=0;
  for(j=0; j<nSpc; j++) {
    if(fuelMoleFrac[j] >= 0.0) {
      fuelSpcIdx[k]=j;
      ++k;
    }
  }
  // record the species indexes found for the oxidizer
  k=0;
  for(j=0; j<nSpc; j++) {
    if(oxidMoleFrac[j] >= 0.0) {
      oxidSpcIdx[k]=j;
      ++k;
    }
  }

  // record the species indexes found for the oxidizer
  k=0;
  for(j=0; j<nSpc; j++) {
    if(traceMoleFrac[j] >= 0.0) {
      traceSpcIdx[k]=j;
      ++k;
    }
  }

  nTrackSpc=0;
  for(j=0; j<nSpc; j++)
    {if(fuelMoleFrac[j] >= 0.0 || oxidMoleFrac[j] >= 0.0 || traceMoleFrac[j] >= 0.0){nTrackSpc++;}}

  //printf("num of tracked species: %d\n",nTrackSpc);
  trackSpcIdx = new int[nTrackSpc];
  nTrackSpc=0;
  for(j=0; j<nSpc; j++)
    {
      if(fuelMoleFrac[j] >= 0.0 || oxidMoleFrac[j] >= 0.0 || traceMoleFrac[j] >= 0.0)
	{trackSpcIdx[nTrackSpc]=j; nTrackSpc++;}
      if(fuelMoleFrac[j] < 0.0){fuelMoleFrac[j]=0.0;}
      if(oxidMoleFrac[j] < 0.0){oxidMoleFrac[j]=0.0;}
      if(traceMoleFrac[j] < 0.0){traceMoleFrac[j]=0.0;}
    }
  normalize(nSpc,fuelMoleFrac);
  normalize(nSpc,oxidMoleFrac);
  // N.B. trace species are NOT normalized
  totalTraceMoleFraction = 0.0;
  for(j=0; j<nSpc; j++) {
    totalTraceMoleFraction += traceMoleFrac[j];
  }

  // at this point the fuel and oxidizer mole fraction compositions should
  // be properly normalized
  fuelOxygenBal=gasMech->getMolarAtomicOxygenRemainder(&fuelMoleFrac[0]);
  oxidOxygenBal=gasMech->getMolarAtomicOxygenRemainder(&oxidMoleFrac[0]);
  moleOxidStoicPerFuel=fabs(fuelOxygenBal)/oxidOxygenBal;

  // declare the remaining mole and mass fraction arrays for setting the
  // initial composition
  freshMoleFrac   = new double[nSpc];
  freshMassFrac   = new double[nSpc];
  exhaustMoleFrac = new double[nSpc];
  exhaustMassFrac = new double[nSpc];
  initMoleFrac    = new double[nSpc];
  initMassFrac    = new double[nSpc];
  // ----- end of parsing the fuel and oxygen mixtures -----
  
  // read simulation constants (non-sweep variables)
  //isConstPres=readFileLine_first_int(inputFptr);
  //if(isConstPres != 0)
  //  {
  //    printf("WARNING: at this time constant pressure IDT is not set up\n");
  //    printf("         calculation will be performed as constant volume\n");
  //    isConstPres=0;
  //  }
  stopTime           = readFileLine_first_double(inputFptr);
  printTime          = readFileLine_first_double(inputFptr);
  maxInternalDt      = readFileLine_first_double(inputFptr);
  maxInternalSteps   = readFileLine_first_int(inputFptr);
  rtol               = readFileLine_first_double(inputFptr);
  atol               = readFileLine_first_double(inputFptr);
  maxCvodeFails1     = readFileLine_first_int(inputFptr);
  maxCvodeFails2     = readFileLine_first_int(inputFptr);
  safetyThresh       = readFileLine_first_double(inputFptr);
  deltaTign          = readFileLine_first_double(inputFptr);
  isContinueAfterIDT = readFileLine_first_int(inputFptr);
  Tref               = readFileLine_first_double(inputFptr);
  
  // read sweep parameters
  nTempRuns  = readFileLine_first_int(inputFptr);
  initTemp = new double[nTempRuns];
  for(j=0; j<nTempRuns; j++)
    {initTemp[j] = readFileLine_first_double(inputFptr);}

  nPresRuns  = readFileLine_first_int(inputFptr);
  initPres = new double[nPresRuns];
  for(j=0; j<nPresRuns; j++)
    {initPres[j] = readFileLine_first_double(inputFptr);}

  nPhiRuns   = readFileLine_first_int(inputFptr);
  phi = new double[nPhiRuns];
  for(j=0; j<nPhiRuns; j++)
    {phi[j] = readFileLine_first_double(inputFptr);}

  nEgrRuns   = readFileLine_first_int(inputFptr);
  egr = new double[nEgrRuns];
  for(j=0; j<nEgrRuns; j++)
    {egr[j] = readFileLine_first_double(inputFptr);}

  nThreshRuns = readFileLine_first_int(inputFptr);
  thresh = new double[nThreshRuns];
  for(j=0; j<nThreshRuns; j++)
    {thresh[j] = readFileLine_first_double(inputFptr);}

  //readFileLine_first_int(inputFptr);
  nKrylovRuns = 1;
  krylovDim = new int[nKrylovRuns];
  for(j=0; j<nKrylovRuns; j++)
    {krylovDim[j] = readFileLine_first_int(inputFptr);}

  //RAW: New options
  doILU = readFileLine_first_int(inputFptr);
  fakeUpdate = readFileLine_first_int(inputFptr);
  threshType = readFileLine_first_int(inputFptr);
  partialPivotThresh = readFileLine_first_double(inputFptr);
  epsLin = readFileLine_first_double(inputFptr);
  nlConvCoeff = readFileLine_first_double(inputFptr);
  oneStepMode = readFileLine_first_int(inputFptr);
  //printAllSteps = readFileLine_first_int(inputFptr);
  permutationType = readFileLine_first_int(inputFptr);

  // set the sweep counters and the initial composition of the first run
  runTotal = nTempRuns*nPresRuns*nPhiRuns*nEgrRuns*nThreshRuns*nKrylovRuns;
  tempId=presId=phiId=egrId=threshId=krylovId=runId=0;
  setInitialComp(phi[phiId],egr[egrId]);

  fclose(inputFptr);
}

idt_sweep_params::~idt_sweep_params()
{
  delete [] krylovDim;
  delete [] thresh;
  delete [] egr;
  delete [] phi;
  delete [] initPres;
  delete [] initTemp;

  delete [] freshMoleFrac;
  delete [] freshMassFrac;
  delete [] exhaustMoleFrac;
  delete [] exhaustMassFrac;
  delete [] initMoleFrac;
  delete [] initMassFrac;

  delete [] trackSpcIdx;
  delete [] fuelSpcIdx;
  delete [] oxidSpcIdx;
  delete [] traceSpcIdx;
  delete [] fuelMoleFrac;
  delete [] oxidMoleFrac;
  delete [] traceMoleFrac;

  delete gasMech;

  free(idtFileName);
  free(logFileName);
  free(thermFileName);
  free(mechFileName);
  free(thistFileName);
}

void idt_sweep_params::setInitialComp(const double phi, const double egr)
{
  int j;

  // compute the fresh mixture mole fraction based on phi and the number
  // of moles of oxidizer needed for stoichiometry per mole of fuel
  for(j=0; j<nSpc; j++)
    {
      freshMoleFrac[j]=phi*fuelMoleFrac[j]
	               +moleOxidStoicPerFuel*oxidMoleFrac[j];
    }
  normalize(nSpc,freshMoleFrac);
  gasMech->getYfromX(freshMoleFrac,freshMassFrac);
  
  if(egr > 0.0)
    {
      // get the ideal exhaust molar composition
      gasMech->getMolarIdealExhaust(freshMoleFrac,exhaustMoleFrac);

      // convert the exhaust mole fractions to mass fractions
      gasMech->getYfromX(exhaustMoleFrac, exhaustMassFrac);
  
      // create the initial mass fraction of the blended intake composition
      // here egr represents the fraction by mass the ideal exhaust composition
      // is in the intake composition
      for(j=0; j<nSpc; j++)
	{initMassFrac[j]=(1.0-egr)*freshMassFrac[j]+egr*exhaustMassFrac[j];}
      // convert the initial intake composition to mole fraction
      gasMech->getXfromY(initMassFrac,initMoleFrac);
    }
  else
    {
      // no egr
      for(j=0; j<nSpc; j++)
	{
	  exhaustMoleFrac[j]=0.0;
	  exhaustMassFrac[j]=0.0;
	  initMoleFrac[j]=freshMoleFrac[j];
	  initMassFrac[j]=freshMassFrac[j];
	}
    }

  //Add in trace species if appropriate
  if(totalTraceMoleFraction > 0.0) {
    for(j=0; j<nSpc; ++j) {
      initMoleFrac[j] *= (1.0-totalTraceMoleFraction);
      initMoleFrac[j] += traceMoleFrac[j];
    }
    normalize(nSpc,initMoleFrac);
    gasMech->getYfromX(initMoleFrac,initMassFrac);
  }
}

void idt_sweep_params::printInitialMoleComp() const
{
  int j;
  printf("# Fuel:\n");
  for(j=0; j<nTrackSpc; j++)
    {
      printf("#      mole of %16s (id = %5d): %.18g\n",
	     gasMech->getSpeciesName(trackSpcIdx[j]),trackSpcIdx[j],
	     fuelMoleFrac[trackSpcIdx[j]]);
    }
  printf("# Oxidizer:\n");
  for(j=0; j<nTrackSpc; j++)
    {
      printf("#      mole of %16s (id = %5d): %.18g\n",
	     gasMech->getSpeciesName(trackSpcIdx[j]),trackSpcIdx[j],
	     oxidMoleFrac[trackSpcIdx[j]]);
    }
  printf("# Trace:\n");
  for(j=0; j<nTrackSpc; j++)
    {
      printf("#      mole of %16s (id = %5d): %.18g\n",
	     gasMech->getSpeciesName(trackSpcIdx[j]),trackSpcIdx[j],
	     traceMoleFrac[trackSpcIdx[j]]);
    }
  printf("# Fresh:\n");
  for(j=0; j<nTrackSpc; j++)
    {
      printf("#      mole of %16s (id = %5d): %.18g\n",
	     gasMech->getSpeciesName(trackSpcIdx[j]),trackSpcIdx[j],
	     freshMoleFrac[trackSpcIdx[j]]);
    }
  printf("# Ideal exhaust:\n");
  for(j=0; j<nTrackSpc; j++)
    {
      printf("#      mole of %16s (id = %5d): %.18g\n",
	     gasMech->getSpeciesName(trackSpcIdx[j]),trackSpcIdx[j],
	     exhaustMoleFrac[trackSpcIdx[j]]);
    }
  printf("# Initial:\n");
  for(j=0; j<nTrackSpc; j++)
    {
      printf("#      mole of %16s (id = %5d): %.18g\n",
	     gasMech->getSpeciesName(trackSpcIdx[j]),trackSpcIdx[j],
	     initMoleFrac[trackSpcIdx[j]]);
    }
}

void idt_sweep_params::incrementRunId()
{
  tempId++;
  if(tempId==nTempRuns)
    {tempId=0; presId++;}
  if(presId==nPresRuns)
    {presId=0; phiId++;}
  if(phiId==nPhiRuns)
    {phiId=0; egrId++;}
  if(egrId==nEgrRuns)
    {egrId=0; threshId++;}
  if(threshId==nThreshRuns)
    {threshId=0; krylovId++;}
  if(krylovId==nKrylovRuns)
    {krylovId=0;} // flip the odometer

  setInitialComp(phi[phiId],egr[egrId]); // update initial composition
  runId++;
}

double idt_sweep_params::getDensity() const
{
  return gasMech->getDensityFromTPY(getInitTemp(),getInitPres(),initMassFrac);
}

void idt_sweep_params::getInitMoleFrac(double x[]) const
{
  int j;
  for(j=0; j<nSpc; j++)
    {x[j]=initMoleFrac[j];}
}

void idt_sweep_params::getInitMassFrac(double y[]) const
{
  int j;
  for(j=0; j<nSpc; j++)
    {y[j]=initMassFrac[j];}
}



double readFileLine_first_double(FILE *fptr)
{
  char *readLine=NULL;
  size_t input_length=0;
  ssize_t read_length;
  int scan_length;
  double val=0.0;
  
  read_length=getline(&readLine,&input_length,fptr);
  if(read_length <= 0)
    {
      printf("WARNING: could not read line in readFileLine_first_double(FILE *)\n");
      printf("         returning 0.0\n");
      fflush(stdout);
      
      if(readLine!=NULL)
	{free(readLine);}

      return 0.0;
    }
  scan_length=sscanf(readLine,"%lf",&val);
  if(scan_length != 1)
    {
      printf("WARNING: could not match data in readFileLine_first_double(FILE *)\n");
      printf("         returning 0.0\n");
      fflush(stdout);
      
      if(readLine!=NULL)
	{free(readLine);}

      return 0.0;
    }

  if(readLine!=NULL)
    {free(readLine);}
 
  return val;
}

int readFileLine_first_int(FILE *fptr)
{
  char *readLine=NULL;
  size_t input_length=0;
  ssize_t read_length;
  int scan_length;
  int val=0;
  
  read_length=getline(&readLine,&input_length,fptr);
  if(read_length <= 0)
    {
      printf("WARNING: could not read line in readFileLine_first_int(FILE *)\n");
      printf("         returning 0\n");
      fflush(stdout);
      
      if(readLine!=NULL)
	{free(readLine);}

      return 0;
    }
  scan_length=sscanf(readLine,"%d",&val);
  if(scan_length != 1)
  {
    printf("WARNING: could not match data in readFileLine_first_int(FILE *)\n");
    printf("         returning 0 (sscanf ret = %d)\n",scan_length);
    fflush(stdout);
  
    if(readLine!=NULL)
  	{free(readLine);}
  
    return 0;
  }

  if(readLine!=NULL)
    {free(readLine);}
 
  return val;
}



char * readFileLine_first_string(FILE *fptr)
{
  const char delimiters[]=" !|\n\t";
  char *readLine=NULL;
  char *token,*str;
  size_t input_length=0;
  ssize_t read_length;
  int scan_length;

  read_length=getline(&readLine,&input_length,fptr);
  if(read_length <= 0)
    {
      printf("WARNING: could not read line in readFileLine_first_string(FILE *, char *)\n");
      printf("         return flag = -1\n");
      fflush(stdout);
      
      if(readLine!=NULL)
	{free(readLine);}

      return NULL;
    } 
  token = strtok(readLine,delimiters);
  scan_length = strlen(token)+1; // add one for terminating byte

  str = (char *)malloc(sizeof(char)*scan_length); // re-allocate string
  strcpy(str,token);

  if(readLine!=NULL)
    {free(readLine);}

 return str;
}

int readFileLine_first_comp(FILE *fptr, zerork::mechanism &mechInp,
			    double moleFrac[])
{
  char specName[MAX_SPEC_NAME_LEN];
  char specMole[MAX_SPEC_MOLE_LEN];
  char *readLine=NULL;
  size_t input_length=0;
  ssize_t read_length;
  int specIdx=-1;
  
  read_length=getline(&readLine,&input_length,fptr);
  if(read_length <= 0)
    {
      printf("WARNING: could not read line in readFileLine_first_comp(FILE *, char *, zerork::mechanism &, double [])\n");
      printf("         return flag = -1\n");
      fflush(stdout);
      
      if(readLine!=NULL)
	{free(readLine);}

      return specIdx;
    }  
  sscanf(readLine,"%s%s",specName,specMole);
  specIdx=mechInp.getIdxFromName(specName);
  moleFrac[specIdx]=atof(specMole);
  //printf("#    found species %s index %d rel. moles %.18g\n",
  //	 specName,specIdx,moleFrac[specIdx]);
  if(readLine!=NULL)
    {free(readLine);}

  return specIdx;
}

double idt_sweep_params::normalize(const int N, double v[])
{
  int j;
  double sum=0.0;
  double invSum=1.0;
  
  for(j=0; j<N; )
    {sum+=v[j]; ++j;}

  if(sum != 0.0)
    {invSum=1.0/sum;}

  for(j=0; j<N; )
    {v[j]*=invSum; ++j;}
  return sum;
}
