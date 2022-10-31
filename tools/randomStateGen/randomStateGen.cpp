#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <utilities/math_utilities.h>

int main(int argc, char *argv[])
{
  int seed;
  int nSpc;
  double Tmin,Tmax,Pmin,Pmax;
  FILE *stateFptr;
  
  int j;

  int isLogT, isLogP;
  double Tstate,Pstate;
  double molFracSum;
  double *molFrac;

  if(argc != 8)
    {
      printf("ERROR: incorrect command line usage\n");
      printf("       use instead %s <seed (-1 = time(0))> <# species>\n",
	     argv[0]);
      printf("                   <min P> <max P> <min T> <max T>\n");
      printf("                   <state file>\n");
      printf("     * Negative P and Ts will use exponential instead of uniform apcing\n");
      exit(-1);
    }
  stateFptr=fopen(argv[7],"w");
  if(stateFptr==NULL)
    {
      printf("ERROR: could not open file %s to write state vector\n",argv[7]);
       exit(-1);
    }
  seed=atoi(argv[1]);
  if(seed==-1)
    {seed=time(0);}
  zerork::utilities::random01seed(seed);

  nSpc=atoi(argv[2]);
  molFrac=(double *)malloc(sizeof(double)*nSpc);
  
  isLogT=isLogP=0;

  Pmin=atof(argv[3]);
  if(Pmin < 0.0)
    {Pmin=-Pmin; isLogP=1;}

  Pmax=atof(argv[4]);
  if(Pmax < 0.0)
    {Pmax=-Pmax; isLogP=1;}

  Tmin=atof(argv[5]);
  if(Tmin < 0.0)
    {Tmin=-Tmin; isLogT=1;}

  Tmax=atof(argv[6]);
  if(Tmax < 0.0)
    {Tmax=-Tmax; isLogT=1;}

  if(isLogP)
    {Pstate=Pmin*exp(log(Pmax/Pmin)*zerork::utilities::random01());}
  else
    {Pstate=Pmin+(Pmax-Pmin)*zerork::utilities::random01();}
      
  if(isLogT)
    {Tstate=Tmin*exp(log(Tmax/Tmin)*zerork::utilities::random01());}
  else
    {Tstate=Tmin+(Tmax-Tmin)*zerork::utilities::random01();}
      
  molFracSum=0.0;
  for(j=0; j<nSpc; j++)
    {
      molFrac[j]=zerork::utilities::random01();
      molFracSum+=molFrac[j];
    }
  molFracSum=1.0/molFracSum;
  for(j=0; j<nSpc; j++)
    {molFrac[j]*=molFracSum;}

  fprintf(stateFptr,"%d\n",nSpc);
  fprintf(stateFptr,"%.18g\n",Pstate);
  fprintf(stateFptr,"%.18g\n",Tstate);
  for(j=0; j<nSpc; j++)
    {fprintf(stateFptr,"%.18g\n",molFrac[j]);}
  fprintf(stateFptr,"# command line:");
  for(j=0; j<argc; j++)
    {fprintf(stateFptr," %s",argv[j]);}
  fprintf(stateFptr,"\n#         seed: %d\n",seed);

  free(molFrac);
  fclose(stateFptr);

  return 0;
}
