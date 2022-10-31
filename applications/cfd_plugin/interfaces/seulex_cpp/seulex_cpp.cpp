

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <assert.h>

#include "seulex_cpp.h"
#include "nvector/nvector_serial.h"
#ifdef ZERORK_CUDA_LIB
#include "nvector/nvector_serial_cuda.h"
#endif

namespace seulex_cpp {

seulex::seulex(int _n, N_Vector _y0)
  :
      n(_n),
      nfcn(0),
      njac(0),
      nstep(0),
      naccept(0),
      nreject(0),
      ndec(0),
      nsol(0),
      nmax(100000),
      km(12),
      nsequ(2),
      lambda(0),
      nrdens(0),
      uround(1.0e-16),
      h(0.0),
      hmax(0.0),
      thet(1.0e-4),
      theta(thet),
      fac1(0.1),
      fac2(4.0),
      fac3(0.7),
      fac4(0.9),
      safe1(0.6),
      safe2(0.93),
      wkfcn(1.0),
      wkjac(5.0),
      wkdec(1.0),
      wksol(1.0),
      wkrow(wkfcn+wksol),
      caljac(false),
      reject(false),
      last(false),
      atov(false),
      err(0.0),
      errold(0.0),
      x(0.0),
      deriv_fcn(NULL),
      jac_fcn(NULL),
      jac_decomp_fcn(NULL),
      jac_solve_fcn(NULL),
      output_fcn(NULL),
      output_fcn_data(NULL),
      user_data(NULL)
{
  //N.B. not storing _y0
  rtol = N_VClone(_y0);
  atol = N_VClone(_y0);
}

seulex::~seulex()
{
  deriv_fcn = NULL;
  jac_fcn = NULL;
  jac_decomp_fcn = NULL;
  jac_solve_fcn = NULL;
  output_fcn = NULL;
  output_fcn_data = NULL;
  user_data = NULL;
  N_VDestroy(rtol);
  N_VDestroy(atol);
}

void seulex::set_nmax(int _nmax)
{
  assert(_nmax >= 0);
  nmax = _nmax;
}

void seulex::set_km(int _km)
{
  assert(_km > 1);
  km = _km;
}

void seulex::set_nsequ(int _nsequ)
{
  assert(_nsequ > 0 && _nsequ < 5);
  nsequ = _nsequ;
}

void seulex::set_lambda(int _lambda)
{
  assert(_lambda == 0 || _lambda == 1);
  lambda=_lambda;
}

void seulex::set_nrdens(int _nrdens)
{
  assert(_nrdens >= 0 && _nrdens <= n);
  nrdens=_nrdens;
}

void seulex::set_uround(double _uround)
{
  assert(_uround > 0. && _uround < 1.);
  uround=_uround;
}

void seulex::set_hinit(double _hinit)
{
  assert(_hinit != 0.);
  h=_hinit;
}

void seulex::set_hmax(double _hmax)
{
  assert(_hmax != 0.);
  hmax=_hmax;
}

void seulex::set_thet(double _thet)
{
  assert(_thet > 0. && _thet < 1.);
  thet=_thet;
}

void seulex::set_step_size_params(double _fac1, double _fac2)
{
  fac1=_fac1;
  fac2=_fac2;
}

void seulex::set_order_params(double _fac3, double _fac4)
{
  fac3=_fac3;
  fac4=_fac4;
}

void seulex::set_safe1safe2(double _safe1, double _safe2)
{
  safe1=_safe1;
  safe2=_safe2;
}

void seulex::set_work_params(double _wkfcn, double _wkjac, double _wkdec, double _wksol)
{
  wkfcn=_wkfcn;
  wkjac=_wkjac;
  wkdec=_wkdec;
  wksol=_wksol;
  wkrow=wkfcn+wksol;
}

void seulex::set_tolerances(N_Vector _rtol, N_Vector _atol)
{
  N_VScale(1.0,_rtol,rtol);
  N_VScale(1.0,_atol,atol);
  check_tolerances();
}

void seulex::set_tolerances(double _rtol, double _atol)
{
  N_VConst(_rtol,rtol);
  N_VConst(_atol,atol);
  check_tolerances();
}

void seulex::check_tolerances()
{
  double min_atol = N_VMin(atol);
  double min_rtol = N_VMin(rtol);
  assert(min_atol > 0.0);
  assert(min_rtol > 10.0*uround);
  thet = std::min(min_rtol,thet);
}


void seulex::set_deriv_fcn(seul_deriv_fcn _deriv_fcn)
{
  deriv_fcn = _deriv_fcn;
}

void seulex::set_jac_fcn(seul_jac_fcn _jac_fcn)
{
  jac_fcn = _jac_fcn;
}

void seulex::set_jac_decomp_fcn(seul_jac_decomp_fcn _jac_decomp_fcn)
{
  jac_decomp_fcn = _jac_decomp_fcn;
}

void seulex::set_jac_solve_fcn(seul_jac_solve_fcn _jac_solve_fcn)
{
  jac_solve_fcn = _jac_solve_fcn;
}

void seulex::set_output_fcn(seul_output_fcn _output_fcn, void* _output_fcn_data)
{
  output_fcn = _output_fcn;
  output_fcn_data = _output_fcn_data;
}

void seulex::set_user_data(void* _user_data)
{
  user_data = _user_data;
}

void seulex::get_integrator_stats(int* _nfcn,int* _njac,int* _nstep,
                            int* _naccept,int* _nreject,int* _ndec,int* _nsol)
{
  *_nfcn = nfcn;
  *_njac = njac;
  *_nstep = nstep;
  *_naccept = naccept;
  *_nreject = nreject;
  *_ndec = ndec;
  *_nsol = nsol;
}

int seulex::solve(double* x_in, double xend, N_Vector y)
{
  //reset counters
  nfcn = 0;
  njac = 0;
  nstep = 0;
  naccept = 0;
  nreject = 0;
  ndec = 0;
  nsol = 0;

  assert(deriv_fcn != NULL);
  assert(jac_fcn != NULL);
  assert(jac_decomp_fcn != NULL);
  assert(jac_solve_fcn != NULL);

  int k, kc, kopt;
  double hopt;

  create_work_vectors(y);
  std::vector<double> hh(km+1,0);
  std::vector<double> w(km+1,0);
  std::vector<double> a(km+1,0);
  std::vector<int> nj(km+1,0);

  setup_step_sequence(nj, a);
  x = *x_in;
  double xstart = x;

  double hmax_save = hmax;
  if(hmax == 0.0) hmax = xend-xstart;
  double posneg = xend-xstart > 0 ? 1.0 : -1.0;
  double min_rtol = N_VMin(rtol);
  double min_atol = N_VMin(atol);
  k = std::max(1,std::min(km-2,int(-log10(min_rtol+min_atol)*0.6+0.5)));
  hmax = std::min(fabs(hmax),fabs(xend-x));
  h = std::max(fabs(h),1.0e-6);
  h = posneg*std::min(h,hmax);
  theta = 2*fabs(thet);
  deriv_fcn(x,y,dy,user_data);
  //if(output_fcn != NULL)
  //{
  //  output_fcn(nstep,x,h,y,dy,output_fcn_data);
  //}

  w[0] = 1.0e30;
  N_VAbs(y,scal);
  N_VProd(rtol,scal,scal);
  N_VLinearSum(1.0, atol, 1.0, scal, scal);
  N_VInv(scal,scal); //N.B. we work with inverse scal.

  caljac=false;
  reject=false;
  //last=false;

g10:
  last=false;
  if(reject) theta=2*fabs(thet);
  atov=false;
  if (0.1*fabs(xend-x)<=fabs(x)*uround) goto g110;
  hopt = h;
  h = posneg*std::min(std::min(fabs(h),fabs(xend-x)),hmax);
  if( (x+1.01*h-xend)*posneg > 0.0 )
  {
    h = xend-x;
    last = true;
  }
  if(theta > thet && !caljac)
  {
    njac += 1;
    deriv_fcn(x,y,dy,user_data);
    jac_fcn(x, y, dy, user_data);
    caljac = true;
  }

//FIRST AND LAST STEP
  if(nstep == 0 || last)
  {
    nstep += 1;
    if (nstep > nmax) goto g120;
    for(int j = 0; j <= k; ++j)
    {
      kc = j;
      seul(x, y, j, hh, w, nj, a);
      if(atov) goto g10;
      if(j > 0 && err <= 1.0) goto g60;
    }
    goto g55;
  }

// --- BASIC INTEGRATION STEP
g30:
  nstep += 1;
  if (nstep > nmax) goto g120;
  kc = k - 1;
  for(int j = 0; j <= kc; ++j)
  {
    seul(x, y, j, hh, w, nj, a);
    if(atov) goto g10;
  }

//C *** *** *** *** *** *** ***
//C --- CONVERGENCE MONITOR
//C *** *** *** *** *** *** ***
  if(k==1 || reject) goto g50;
  if(err <= 1.0) goto g60;
  if(err > ((double) (nj[k+1]*nj[k]))*4.) goto g100;

g50:
  seul(x, y, k, hh, w, nj, a);
  if(atov) goto g10;
  kc = k;
  if(err <= 1.0) goto g60;

g55:
  if(err > ((double) nj[k+1])*2.) goto g100;
  kc = k + 1;
// --- HOPE FOR CONVERGENCE IN LINE K+1
  seul(x, y, kc, hh, w, nj, a);
  if(atov) goto g10;

  if(err > 1.0) goto g100;

//C*** *** *** *** *** *** ***
//C --- STEP IS ACCEPTED
//C *** *** *** *** *** *** ***
g60:
  x = x + h;
  //N.B. : could be improved by fusing (also above)
  N_VAbs(t[0],scal);
  N_VProd(rtol,scal,scal);
  N_VLinearSum(1.0, atol, 1.0, scal, scal);
  N_VInv(scal,scal);
  N_VScale(1.0,t[0],y);
  naccept += 1;
  caljac = false;
  //TODO: Interpolation functions for continous output
  if(output_fcn != NULL)
  {
    int ofcn_flag = output_fcn(naccept,x,h,y,dy,output_fcn_data);
    if(ofcn_flag != 0) {
      goto g110;
    }
  }
//C --- COMPUTE OPTIMAL ORDER
  if(kc == 1)
  {
    kopt = std::min(2,km-1);
    if (reject) kopt=1;
    goto g80;
  }
  if(kc <= k)
  {
    kopt=kc;
    if (w[kc-1] < w[kc]*fac3) kopt=kc-1;
    if (w[kc] < w[kc-1]*fac4) kopt=std::min(kc+1,km-1);
  }
  else
  {
    kopt=kc-1;
    if (kc > 2 && w[kc-2] < w[kc-1]*fac3) kopt=kc-2;
    if (w[kc] < w[kopt]*fac4) kopt=std::min(kc,km-1);
  }
  //printf("x, h, k, kopt, kc = %g, %g, %d, %d, %d\n",x,h,k,kopt,kc);

g80:
  if (reject)
  {
     k=std::min(kopt,kc);
     h=posneg*std::min(fabs(h),fabs(hh[k]));
     reject=false;
     goto g10;
  }
//C --- COMPUTE STEP SIZE FOR NEXT STEP
  if (kopt <= kc)
  {
    h=hh[kopt];
  }
  else
  {
    if (kc < k && w[kc] < w[kc-1]*fac4)
    {
      h=hh[kc]*a[kopt+1]/a[kc];
    }
    else
    {
      h=hh[kc]*a[kopt]/a[kc];
    }
  }
  k=kopt;
  h=posneg*fabs(h);
  goto g10;

//C *** *** *** *** *** *** ***
//C --- STEP IS REJECTED
//C *** *** *** *** *** *** ***
g100:
  k=std::min(k,std::min(kc,km-1));
  if (k > 2 && w[k-1] < w[k]*fac3) k=k-1;
  nreject += 1;
  h=posneg*hh[k];
  last=false;
  reject=true;
  if (caljac) goto g30;
  goto g10;

//C --- SOLUTION EXIT
g110:
  h=hopt;
  destroy_work_vectors();
  hmax = hmax_save;
  *x_in = x;
  return 0;

//C --- FAIL EXIT
g120:
  printf("  EXIT OF SEULEX AT X=%g   H=%g  nstep=%d\n",x,h,nstep);
  destroy_work_vectors();
  hmax = hmax_save;
  *x_in = x;
  return 1;
}

//C *** *** *** *** *** *** ***
//C     S U B R O U T I N E    S E U L
//C *** *** *** *** *** *** ***
int seulex::seul(double x,
             N_Vector y,
             int jj,
             std::vector<double>& hh,
             std::vector<double>& w,
             std::vector<int>& nj,
             std::vector<double>& a)
{
  int m;
  double fac,facmin,expo;
//C --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE
//C --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATE
//C --- OF THE OPTIMAL STEP SIZE
  double hj=h/nj[jj];
  double hji=1.0/hj;
//C *** *** *** *** *** *** ***
//C  COMPUTE THE MATRIX E AND ITS DECOMPOSITION
//C *** *** *** *** *** *** ***
  int ier = jac_decomp_fcn(hj,user_data);
  ndec += 1;
  if(ier != 0)
  {
//    printf("jac_decomp_fcn ier=%d\n",ier);
    goto g79;
  }
//C *** *** *** *** *** *** ***
//C --- STARTING PROCEDURE
//C *** *** *** *** *** *** ***
  ier = deriv_fcn(x+hj,y,dy,user_data);
  nfcn += 1;
  if(ier != 0) goto g79;
  N_VScale(1.0,y,yh);
  N_VScale(1.0,dy,del);
  //printf("nstep, nreject, njac, ndec, nsol, jj, hji, thet, theta = %d, %d, %d, %d, %d, %d, %g, %g, %g\n",
  //       nstep,nreject,njac,ndec,nsol,jj,hji,thet,theta);
  //DEBUG N_VPrint_Serial(del);
  jac_solve_fcn(x, y, dy, del, tmp1, user_data);
  N_VScale(1.0, tmp1, del);
  N_VScale(hj,del,del);

  nsol += 1;
  m=nj[jj];


//C *** *** *** *** *** *** ***
//C --- SEMI-IMPLICIT EULER METHOD
//C *** *** *** *** *** *** ***
  if (m > 1)
  {
    for(int mm = 1; mm < m; ++mm)
    {
      N_VLinearSum(1.0,yh,1.0,del,yh);

      ier = deriv_fcn(x+hj*(mm+1),yh,dyh,user_data);
      nfcn+=1;
      if(ier != 0) goto g79;
      if (mm == 1 && jj <= 1) //mm 1-based; jj 0-based
      {
//C --- STABILITY CHECK
         double del1 = N_VWL2Norm(del,scal);
         if(isnan(del1)) goto g79;
         ier = deriv_fcn(x+hj,yh,wh,user_data);
         nfcn+=1;
         if(ier != 0) goto g79;
         N_VLinearSum(1.0,wh,-hji,del,del);
         jac_solve_fcn(x+hj, yh, wh, del, tmp1, user_data);
         N_VScale(1.0, tmp1, del);
         N_VScale(hj,del,del);
         nsol+=1;

         double max_del = N_VMaxNorm(del);
         if(max_del > 1e15) goto g79;

         double del2 = N_VWL2Norm(del,scal);
         if(isnan(del2)) goto g79;
         theta=del2/std::max(1.0,del1);
         //theta=del2/std::min(1.0,del1+1.0e-30);
         if (theta > 1.0) goto g79;
      }
      //N_VScale(hj,dyh,dyh);
      jac_solve_fcn(x+hj*(mm+1), yh, dyh, dyh, tmp1, user_data);
      N_VScale(1.0, tmp1, dyh);
      N_VScale(hj,dyh,dyh);
      nsol+=1;
      N_VScale(1.0,dyh,del);
//      IF (IOUT.EQ.2.AND.MM.GE.M-JJ) THEN
//        IPT=IPT+1
//        DO I=1,NRD
//          FSAFE(IPT,I)=DEL(ICOMP(I))
//        END DO
//      END IF
    }
  }
  N_VLinearSum(1.0,yh,1.0,del,t[jj]);
  //DEBUG N_VPrint_Serial(t[jj]);
//C *** *** *** *** *** *** ***
//C --- POLYNOMIAL EXTRAPOLATION
//C *** *** *** *** *** *** ***
  if(jj==0) return 0;
  for(int l = jj; l > 0; --l)
  {
     //N.B. OpenFOAM implementation uses a table
     double fac=((double) nj[jj])/((double)nj[l-1])-1.0;
     fac = 1.0/fac;
     N_VLinearSum(1.0,t[l],-1.0,t[l-1],t[l-1]);
     N_VLinearSum(1.0,t[l],fac,t[l-1],t[l-1]);
  }
  err = 0.0;
  N_VLinearSum(1.0, t[0], -1.0, t[1], wh); //USING wh for temp storage here
  N_VAbs(wh,wh);
  err = N_VWrmsNorm(scal,wh);
  if(isnan(err)) goto g79;
  if(err < 0.0) goto g79;
  if (err > 1.0e15) goto g79;
  if (jj > 1 && err >= errold) goto g79;
  errold = std::max(4*err,1.0);
//C --- COMPUTE OPTIMAL STEP SIZES
  expo=1.0/(jj+1);
  facmin=pow(fac1,expo);
  fac=std::min(fac2/facmin,std::max(facmin,pow((err/safe1),expo)/safe2));
  fac=1.0/fac;
  hh[jj]=std::min(fabs(h)*fac,hmax);
  w[jj]=a[jj]/hh[jj];
  return 0;
g79:
  atov=true;
  h=h*0.5;
  reject=true;
  return 1;
};

void seulex::setup_step_sequence(std::vector<int>& nj, std::vector<double>& a)
{
// --- DEFINE THE STEP SIZE SEQUENCE
  if (nsequ == 1)
  {
    nj[0]=1;
    nj[1]=2;
    nj[2]=3;
    for(int i = 3; i <= km; ++i)
    {
      nj[i]=2*nj[i-2];
    }
  }
  else if (nsequ == 2)
  {
    nj[0]=2;
    nj[1]=3;
    for(int i = 2; i <= km; ++i)
    {
      nj[i]=2*nj[i-2];
    }
  }
  else if (nsequ == 3)
  {
    for(int i = 0; i <= km; ++i)
    {
      nj[i]=i+1;
    }
  }
  else if (nsequ == 4)
  {
    for(int i = 0; i <= km; ++i)
    {
      nj[i]=i+2;
    }
  }
  a[0] = wkjac+nj[0]*wkrow+wkdec;
  for(int i = 1; i <= km; ++i)
  {
    a[i]=a[i-1]+(nj[i]-1)*wkrow+wkdec;
  }
}

void seulex::create_work_vectors(N_Vector y)
{
  yh = N_VClone(y);
  dy = N_VClone(y);
  dyh = N_VClone(y);
  del = N_VClone(y);
  wh = N_VClone(y);
  scal = N_VClone(y);
  tmp1 = N_VClone(y);
  t.resize(km+1);
  for(int i = 0; i <= km; ++i)
  {
    t[i] = N_VClone(y);
  }
}

void seulex::destroy_work_vectors()
{
  N_VDestroy(yh);
  N_VDestroy(dy);
  N_VDestroy(dyh);
  N_VDestroy(del);
  N_VDestroy(wh);
  N_VDestroy(scal);
  N_VDestroy(tmp1);
  for(int i = 0; i <= km; ++i)
  {
    N_VDestroy(t[i]);
  }
}

} //end namespace seulex_cpp
