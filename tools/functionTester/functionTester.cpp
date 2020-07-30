#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "zerork/mechanism.h"
#include "functionTester.h"

int main(int argc, char *argv[])
{
  FILE *stateFptr,*outputFptr;
  zerork::mechanism *mech;

  int j;

  int nSpc, nRxn, nmechSpc,nStep;
  char readLine[MAX_LINE_LEN];
  double *moleFrac;
  double *massFrac,*conc,*cvSpc,*cpSpc,*uSpc,*hSpc,*gSpc,*molwtSpc;
  double *Kfwd,*Krev,*stepROP;
  double *netSpc,*createSpc,*destroySpc;

  double *netSpcInternal;

  double pres,Temp,rvol,dens,molwtMix,presConvert,cpMix,cvMix,uMix,hMix;
  double gasConstant;

  if(argc != 6)
    {
      printf("ERROR: incorrect command line usage.\n");
      printf("       use instead %s <ck2 mech file>  <ck2 thermo file> <ck2 converter output file> <state file> <output file>\n",argv[0]);
      exit(-1);
    }

  stateFptr=fopen(argv[4],"r");
  if(stateFptr==NULL)
    {
      printf("ERROR: could not open state vector file %s for read\n",
	     argv[4]);
      exit(-1);
    }
  outputFptr=fopen(argv[5],"w");
  if(outputFptr==NULL)
    {
      printf("ERROR: could not open output file %s for write\n",
	     argv[5]);
      exit(-1);
    }
  mech=new zerork::mechanism(argv[1],argv[2],argv[3]);

  nmechSpc=mech->getNumSpecies();
  nRxn=mech->getNumReactions();
  nStep=mech->getNumSteps();
  // parse the input state vector
  GetLine(stateFptr,readLine,'\n',MAX_LINE_LEN);
  sscanf(readLine,"%d",&nSpc);
  if(nSpc != nmechSpc)
    {
      printf("ERROR: number of species in mechanism file %d\n",nmechSpc);
      printf("       differs from state file %d\n",nSpc);
      exit(-1);
    }
  
  GetLine(stateFptr,readLine,'\n',MAX_LINE_LEN);
  sscanf(readLine,"%lf",&pres);
  GetLine(stateFptr,readLine,'\n',MAX_LINE_LEN);
  sscanf(readLine,"%lf",&Temp);

  moleFrac   = new double[nSpc];
  massFrac   = new double[nSpc];
  conc       = new double[nSpc];
  cvSpc      = new double[nSpc];
  cpSpc      = new double[nSpc];
  uSpc       = new double[nSpc];
  hSpc       = new double[nSpc];
  gSpc       = new double[nSpc];
  molwtSpc   = new double[nSpc];
  Kfwd       = new double[nRxn];
  Krev       = new double[nRxn];
  stepROP    = new double[nStep];
  netSpc     = new double[nSpc];
  netSpcInternal = new double[nSpc];
  createSpc  = new double[nSpc];
  destroySpc = new double[nSpc];

  for(j=0; j<nSpc-1; j++)
    {
      GetLine(stateFptr,readLine,'\n',MAX_LINE_LEN);
      sscanf(readLine,"%lf",&moleFrac[j]);
    }
  // readin the last mole frac ignoring the carriage return
  fscanf(stateFptr,"%lf",&moleFrac[nSpc-1]);
 
  // first echo the input
  fprintf(outputFptr,"%24d ! [#]             number of species\n",nSpc);
  fprintf(outputFptr,"%32.24e ! [Pa]            input pressure\n",pres);
  fprintf(outputFptr,"%32.24e ! [K]             input temperature\n",Temp);
  for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [-]             input mole frac - %s\n",
	      moleFrac[j],mech->getSpeciesName(j));
    }
  // compute the mass fractions
  mech->getYfromX(moleFrac,massFrac);

  // compute the mixture properties
  dens=mech->getDensityFromTPY(Temp,pres,massFrac);
  rvol=1.0/dens;
  molwtMix=mech->getMolWtMixFromY(massFrac);
  presConvert=mech->getPressureFromTVY(Temp,rvol,massFrac);
  uMix=mech->getMassIntEnergyFromTY(Temp,massFrac,uSpc);
  hMix=mech->getMassEnthalpyFromTY(Temp,massFrac,hSpc);
  cpMix=mech->getMassCpFromTY(Temp,massFrac,cpSpc);
  cvMix=mech->getMassCvFromTY(Temp,massFrac,cvSpc);
  gasConstant=mech->getGasConstant();

  fprintf(outputFptr,"%32.24e ! [kg/m^3]        density\n",dens);
  fprintf(outputFptr,"%32.24e ! [m^3/kg]        relative volume\n",rvol);
  fprintf(outputFptr,"%32.24e ! [kg/kmol]       molecular weight\n",molwtMix);
  fprintf(outputFptr,"%32.24e ! [Pa]            pressure (from rel volume)\n",
	  presConvert);
  fprintf(outputFptr,"%32.24e ! [J/kg]          mixture internal energy\n",uMix);
  fprintf(outputFptr,"%32.24e ! [J/kg]          mixture enthalpy\n",hMix);
  fprintf(outputFptr,"%32.24e ! [J/kg/K]        specific heat (const vol)\n",
	  cvMix);
  fprintf(outputFptr,"%32.24e ! [J/kg/K]        specific heat (const pres)\n",
	  cpMix);
  fprintf(outputFptr,"%32.24e ! [J/kmol/K]      univerisal gas const\n",
	  gasConstant);

  // calculate species properties
  mech->getMolWtSpc(molwtSpc);
  mech->getCfromVY(rvol,massFrac,conc);
  mech->getNonDimGibbsFromT(Temp,gSpc);

  // calculate the reaction rates
  mech->getReactionRates(Temp,conc,netSpc,createSpc,destroySpc,stepROP);

  // calculate the net reaction rate using the internal reation rate workspace  
  mech->getNetReactionRates(Temp,conc,netSpcInternal);

  // perform an internal check of the two methods for the net species rates
  double maxDiff=0.0;
  for(j=0; j<nSpc; j++) {

    if(fabs(netSpcInternal[j]-netSpc[j]) > maxDiff) {
      maxDiff=fabs(netSpcInternal[j]-netSpc[j]);
    }
  }
  if(maxDiff == 0.0) {
    printf("# internal check for net species production rate: passed\n");
  }
  else {
    printf("# internal check for net species production rate: failed\n");
    printf("# species id: getReactionRates(...)  getNetReactionRates(...)  |delta|\n");
    for(j=0; j<nSpc; j++) {
      printf("%d: %.18g  %.18g  %.18g\n",
             j,
             netSpc[j],
             netSpcInternal[j],
             fabs(netSpcInternal[j]-netSpc[j]));
    }
  }

  //for(j=0; j<nSpc; j++)
  // {netSpc[j]=createSpc[j]=destroySpc[j]=0.0;}
  
  // calculate reaction rate constants
  mech->getKrxnFromTC(Temp,conc,Kfwd,Krev);

  for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [kg/kmol]       molecular weight - %s\n",
	      molwtSpc[j],mech->getSpeciesName(j));
    }
  for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [-]             mass frac - %s\n",
	      massFrac[j],mech->getSpeciesName(j));
    }
  for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [kmol/m^3]      concentration - %s\n",
	      conc[j],mech->getSpeciesName(j));
    }
 
 for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [J/kg]          internal energy - %s\n",
	      uSpc[j],mech->getSpeciesName(j));
    }  

 for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [J/kg]          enthalpy - %s\n",
	      hSpc[j],mech->getSpeciesName(j));
    }  

 for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [J/kg/K]        specific heat (const vol) - %s\n",cvSpc[j],mech->getSpeciesName(j));
    }  

 for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [J/kg/K]        specific heat (const pres) - %s\n",cpSpc[j],mech->getSpeciesName(j));
    }  

 for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [-]             (S/R - H/RT) - %s\n",-gSpc[j],mech->getSpeciesName(j));
    } 
 for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [kmol/m^3/s]     net production rate - %s\n",netSpc[j],mech->getSpeciesName(j));
    }
 for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [kmol/m^3/s]     creation rate - %s\n",createSpc[j],mech->getSpeciesName(j));
    }
 for(j=0; j<nSpc; j++)
    {
      fprintf(outputFptr,"%32.24e ! [kmol/m^3/s]     destruction rate - %s\n",destroySpc[j],mech->getSpeciesName(j));
    }
 for(j=0; j<nRxn; j++)
   {
      fprintf(outputFptr,"%32.24e ! [(m^3/kmol)^n/s] fwd rate constant rxn - %6d\n",
	     Kfwd[j],j+1);
   }
 for(j=0; j<nRxn; j++)
   {
     fprintf(outputFptr,"%32.24e ! [(m^3/kmol)^n/s] rev rate constant rxn - %6d\n",
	     Krev[j],j+1);
   }

  delete mech;
  delete [] moleFrac;
  delete [] massFrac;
  delete [] conc;
  delete [] cvSpc;
  delete [] cpSpc;
  delete [] uSpc;
  delete [] hSpc;
  delete [] gSpc;
  delete [] molwtSpc;
  delete [] stepROP;
  delete [] netSpc;
  delete [] netSpcInternal;
  delete [] createSpc;
  delete [] destroySpc;
  delete [] Kfwd;
  delete [] Krev;

  fclose(outputFptr);
  fclose(stateFptr);

  return 0;
}

int GetLine(FILE *InFile,char *ReadLine, char UntilChar, int MaxChar)
{
  int CharCtr=0;
  char LocalReadChar='\0';


  while((LocalReadChar != UntilChar) && CharCtr < (MaxChar-1))
    {
      fscanf(InFile,"%c",&LocalReadChar);
      //printf("LocalReadChar[%d]: %c    UntilChar: %c\n",CharCtr,
      //             LocalReadChar,UntilChar);
      ReadLine[CharCtr]=LocalReadChar;
      CharCtr++;
    }
  if(CharCtr == (MaxChar-1) && LocalReadChar != UntilChar) // ran out of space
    {ReadLine[0]='\0'; return -1;}

  ReadLine[CharCtr]='\0';
  return CharCtr; // exit normally
}
