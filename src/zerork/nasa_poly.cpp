#include "nasa_poly.h"
#include "constants.h"
#include <stdio.h>
#include <math.h>
#include <vector>
#include <assert.h>

namespace zerork {

nasa_poly_group::nasa_poly_group(const int inpSpc,
  double const * const inpCoef, double const * const inpTlow,
  double const * const inpThigh)
{
  int j,k;

//  if(inpSpc < 1) {
//    printf("ERROR: allocator nasa_poly_group(...) failed\n");
//    printf("       because the number of species in the group %d < 1.\n",
//           inpSpc);
//  }
  assert(inpSpc > 0);
  thermoCoef = new double[inpSpc*LDA_THERMO_POLY_D5R2];
  Tlow = new double[inpSpc];
  Thigh = new double[inpSpc];
  // MJM future fix - add allocation checks

  // intialize the thermo group data
  nGroupSpc=inpSpc;
  for(j=0; j<inpSpc; j++)
  {
    for(k=0; k<NUM_THERMO_POLY_D5R2; k++)
    {
      thermoCoef[j*LDA_THERMO_POLY_D5R2+k]=
        inpCoef[j*NUM_THERMO_POLY_D5R2+k];
    }
    for(k=NUM_THERMO_POLY_D5R2; k<LDA_THERMO_POLY_D5R2; k++)
    {thermoCoef[j*LDA_THERMO_POLY_D5R2+k]=0.0;}
    Tlow[j]=inpTlow[j];
    Thigh[j]=inpThigh[j];
  }
}

nasa_poly_group::~nasa_poly_group()
{
  delete [] thermoCoef;
  delete [] Tlow;
  delete [] Thigh;
}

void nasa_poly_group::getCp_R(const double T, double Cp_R[]) const
{
  int j,coefAddr;
  double Tmid;

  coefAddr=0;
  for(j=0; j<nGroupSpc; j++)
    {
      Tmid=thermoCoef[coefAddr];
      if(T < Tmid)
	{
	  Cp_R[j]=     thermoCoef[coefAddr+1]+
	            T*(thermoCoef[coefAddr+2]+
		    T*(thermoCoef[coefAddr+3]+
		    T*(thermoCoef[coefAddr+4]+
		    T* thermoCoef[coefAddr+5])));
	}
      else
	{
	  Cp_R[j]=     thermoCoef[coefAddr+8 ]+
	            T*(thermoCoef[coefAddr+9 ]+
		    T*(thermoCoef[coefAddr+10]+
		    T*(thermoCoef[coefAddr+11]+
		    T* thermoCoef[coefAddr+12])));
	}
      coefAddr+=LDA_THERMO_POLY_D5R2;
    }
}

void nasa_poly_group::getH_RT(const double T, double H_RT[]) const
{
  int j,coefAddr;
  double Tmid;
  double invT=1.0/T;
  double hMult[4];

  coefAddr=0;

  hMult[0]=T;
  hMult[1]=T*hMult[0];
  hMult[2]=T*hMult[1];
  hMult[3]=T*hMult[2];

  hMult[0]*=0.50000000000000000000; // h0 = T/2
  hMult[1]*=0.33333333333333333333; // h1 = T^2/3
  hMult[2]*=0.25000000000000000000; // h2 = T^3/4
  hMult[3]*=0.20000000000000000000; // h3 = T^4/5

  for(j=0; j<nGroupSpc; j++)
  {
    Tmid=thermoCoef[coefAddr];
    if(T < Tmid)
    {
      H_RT[j]=thermoCoef[coefAddr+1]+
        thermoCoef[coefAddr+2]*hMult[0]+
        thermoCoef[coefAddr+3]*hMult[1]+
        thermoCoef[coefAddr+4]*hMult[2]+
        thermoCoef[coefAddr+5]*hMult[3]+
        thermoCoef[coefAddr+6]*invT;
    }
    else
    {
      H_RT[j]=thermoCoef[coefAddr+8 ]+
        thermoCoef[coefAddr+9 ]*hMult[0]+
        thermoCoef[coefAddr+10]*hMult[1]+
        thermoCoef[coefAddr+11]*hMult[2]+
        thermoCoef[coefAddr+12]*hMult[3]+
        thermoCoef[coefAddr+13]*invT;
    }
    coefAddr+=LDA_THERMO_POLY_D5R2;
  }
}
//
// calculate the non-dimensional Gibbs free energy at the standard state
// pressure 1 atm
//
// G/RT = H/RT - S/R
void nasa_poly_group::getG_RT(const double T, double G_RT[]) const
{
  int j,coefAddr;
  double Tmid;
  double invT=1.0/T;
  double gMult[5];

  gMult[0]=1.0-log(T);
  gMult[1]=T;
  gMult[2]=T*gMult[1];
  gMult[3]=T*gMult[2];
  gMult[4]=T*gMult[3];

                                      // g0 = 1 - ln(T)
  gMult[1]*=-0.50000000000000000000;  // g1 = - T/2
  gMult[2]*=-0.16666666666666666667;  // g2 = - T^2/6
  gMult[3]*=-0.08333333333333333333;  // g3 = - T^3/12
  gMult[4]*=-0.05000000000000000000;  // g4 = - T^4/20

  coefAddr=0;
  for(j=0; j<nGroupSpc; j++)
    {
      Tmid=thermoCoef[coefAddr];
      if(T < Tmid)
	{
	  G_RT[j]=thermoCoef[coefAddr+1]*gMult[0]+
	          thermoCoef[coefAddr+2]*gMult[1]+
	          thermoCoef[coefAddr+3]*gMult[2]+
	          thermoCoef[coefAddr+4]*gMult[3]+
	          thermoCoef[coefAddr+5]*gMult[4]+
	          thermoCoef[coefAddr+6]*invT-
	          thermoCoef[coefAddr+7];
	}
      else
	{
	  G_RT[j]=thermoCoef[coefAddr+8 ]*gMult[0]+
	          thermoCoef[coefAddr+9 ]*gMult[1]+
	          thermoCoef[coefAddr+10]*gMult[2]+
	          thermoCoef[coefAddr+11]*gMult[3]+
	          thermoCoef[coefAddr+12]*gMult[4]+
	          thermoCoef[coefAddr+13]*invT-
	          thermoCoef[coefAddr+14];
	}
      coefAddr+=LDA_THERMO_POLY_D5R2;
    }
}


void nasa_poly_group::getCp_R_mr(const int nReactors, const double T[], double Cp_R[]) const
{
  int j,k,coefAddr;
  double Tmid;

  coefAddr=0;
  for(j=0; j<nGroupSpc; j++)
    {
      Tmid=thermoCoef[coefAddr];
      for(k=0;k<nReactors;++k)
      {
          if(T[k] < Tmid)
	    {
	      Cp_R[nReactors*j+k]=     thermoCoef[coefAddr+1]+
	                T[k]*(thermoCoef[coefAddr+2]+
		        T[k]*(thermoCoef[coefAddr+3]+
		        T[k]*(thermoCoef[coefAddr+4]+
		        T[k]* thermoCoef[coefAddr+5])));
	    }
          else
	    {
	      Cp_R[nReactors*j+k]=     thermoCoef[coefAddr+8 ]+
	                T[k]*(thermoCoef[coefAddr+9 ]+
		        T[k]*(thermoCoef[coefAddr+10]+
		        T[k]*(thermoCoef[coefAddr+11]+
		        T[k]* thermoCoef[coefAddr+12])));
	    }
      }
      coefAddr+=LDA_THERMO_POLY_D5R2;
    }
}

void nasa_poly_group::getH_RT_mr(const int nReactors, const double T[], double H_RT[]) const
{
  int j,k,coefAddr;
  double Tmid;
  std::vector<double> invT(nReactors);
  double hMult[4];

  for(k=0;k<nReactors;++k)
  {
      invT[k] =1.0/T[k];
  }
  coefAddr=0;

  for(j=0; j<nGroupSpc; j++)
    {
      Tmid=thermoCoef[coefAddr];
      for(k=0;k<nReactors;++k)
      {
          hMult[0]=T[k];
          hMult[1]=T[k]*hMult[0];
          hMult[2]=T[k]*hMult[1];
          hMult[3]=T[k]*hMult[2];

          hMult[0]*=0.50000000000000000000; // h0 = T/2
          hMult[1]*=0.33333333333333333333; // h1 = T^2/3
          hMult[2]*=0.25000000000000000000; // h2 = T^3/4
          hMult[3]*=0.20000000000000000000; // h3 = T^4/5
          if(T[k] < Tmid)
	    {
	      H_RT[nReactors*j+k]=thermoCoef[coefAddr+1]+
	              thermoCoef[coefAddr+2]*hMult[0]+
		      thermoCoef[coefAddr+3]*hMult[1]+
		      thermoCoef[coefAddr+4]*hMult[2]+
		      thermoCoef[coefAddr+5]*hMult[3]+
                      thermoCoef[coefAddr+6]*invT[k];
	    }
          else
	    {
	      H_RT[nReactors*j+k]=thermoCoef[coefAddr+8 ]+
	              thermoCoef[coefAddr+9 ]*hMult[0]+
		      thermoCoef[coefAddr+10]*hMult[1]+
		      thermoCoef[coefAddr+11]*hMult[2]+
		      thermoCoef[coefAddr+12]*hMult[3]+
                      thermoCoef[coefAddr+13]*invT[k];
	    }
      }
      coefAddr+=LDA_THERMO_POLY_D5R2;
    }
}

void nasa_poly_group::getThermoCoeffs(double coeffs[]) const
{
  int i,j;
  for(j=0; j<nGroupSpc; j++) {
    for(i=0; i<LDA_THERMO_POLY_D5R2; i++)
    {
      coeffs[j*LDA_THERMO_POLY_D5R2 + i] = thermoCoef[j*LDA_THERMO_POLY_D5R2 + i];
    }
  }
}

} // namespace zerork
