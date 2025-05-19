

#include "sodex_cpp.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <algorithm>

namespace sodex_cpp {

sodex::sodex(int _n, N_Vector _y0)
    : n(_n),
      nfcn(0),
      nstep(0),
      naccept(0),
      nreject(0),
      ndec(0),
      nsol(0),
      nmax(100000),
      km(6),
      nsequ(1),
      nrdens(0),
      uround(1.0e-16),  // TODO:
      h(0.0),
      hmax(0.0),
      hmaxn(0.0),
      thet(1.0e-4),
      theta(thet),
      fac1(0.1),
      fac2(4.0),
      fac3(0.9),
      fac4(0.9),
      safe1(0.8),
      safe2(0.93),
      wkfcn(1.0),
      wkjac(5.0),
      wkdec(1.0),
      wksol(1.0),
      wkrow(wkfcn + wksol),
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
      user_data(NULL) {
  // N.B. not storing _y0
  rtol = N_VClone(_y0);
  atol = N_VClone(_y0);
}

sodex::~sodex() {
  deriv_fcn = NULL;
  jac_fcn = NULL;
  jac_decomp_fcn = NULL;
  jac_solve_fcn = NULL;
  output_fcn = NULL;
  user_data = NULL;
  output_fcn_data = NULL;
  N_VDestroy(rtol);
  N_VDestroy(atol);
}

void sodex::set_nmax(int _nmax) {
  assert(_nmax >= 0);
  nmax = _nmax;
}

void sodex::set_km(int _km) {
  assert(_km > 1);
  km = _km;
}

void sodex::set_nsequ(int _nsequ) {
  assert(_nsequ > 0 && _nsequ < 2);
  nsequ = _nsequ;
}

void sodex::set_nrdens(int _nrdens) {
  assert(_nrdens >= 0 && _nrdens <= n);
  nrdens = _nrdens;
}

void sodex::set_uround(double _uround) {
  assert(_uround > 0. && _uround < 1.);
  uround = _uround;
}

void sodex::set_hinit(double _hinit) {
  assert(_hinit != 0.);
  h = _hinit;
}

void sodex::set_hmax(double _hmax) {
  assert(_hmax != 0.);
  hmax = _hmax;
}

void sodex::set_thet(double _thet) {
  assert(_thet > 0. && _thet < 1.);
  thet = _thet;
}

void sodex::set_step_size_params(double _fac1, double _fac2) {
  fac1 = _fac1;
  fac2 = _fac2;
}

void sodex::set_order_params(double _fac3, double _fac4) {
  fac3 = _fac3;
  fac4 = _fac4;
}

void sodex::set_safe1safe2(double _safe1, double _safe2) {
  safe1 = _safe1;
  safe2 = _safe2;
}

void sodex::set_work_params(double _wkfcn, double _wkjac, double _wkdec,
                            double _wksol) {
  wkfcn = _wkfcn;
  wkjac = _wkjac;
  wkdec = _wkdec;
  wksol = _wksol;
  wkrow = wkfcn + wksol;
}

void sodex::set_tolerances(N_Vector _rtol, N_Vector _atol) {
  N_VScale(1.0, _rtol, rtol);
  N_VScale(1.0, _atol, atol);
  check_tolerances();
}

void sodex::set_tolerances(double _rtol, double _atol) {
  N_VConst(_rtol, rtol);
  N_VConst(_atol, atol);
  check_tolerances();
}

void sodex::check_tolerances() {
  double min_atol = N_VMin(atol);
  double min_rtol = N_VMin(rtol);
  assert(min_atol > 0.0);
  assert(min_rtol > 10.0 * uround);
}

void sodex::set_deriv_fcn(sodex_deriv_fcn _deriv_fcn) {
  deriv_fcn = _deriv_fcn;
}

void sodex::set_jac_fcn(sodex_jac_fcn _jac_fcn) { jac_fcn = _jac_fcn; }

void sodex::set_jac_decomp_fcn(sodex_jac_decomp_fcn _jac_decomp_fcn) {
  jac_decomp_fcn = _jac_decomp_fcn;
}

void sodex::set_jac_solve_fcn(sodex_jac_solve_fcn _jac_solve_fcn) {
  jac_solve_fcn = _jac_solve_fcn;
}

void sodex::set_output_fcn(sodex_output_fcn _output_fcn,
                           void* _output_fcn_data) {
  output_fcn = _output_fcn;
  output_fcn_data = _output_fcn_data;
}

void sodex::set_user_data(void* _user_data) { user_data = _user_data; }

void sodex::get_integrator_stats(int* _nfcn, int* _njac, int* _nstep,
                                 int* _naccept, int* _nreject, int* _ndec,
                                 int* _nsol) {
  *_nfcn = nfcn;
  *_njac = njac;
  *_nstep = nstep;
  *_naccept = naccept;
  *_nreject = nreject;
  *_ndec = ndec;
  *_nsol = nsol;
}

int sodex::solve(double* x_in, double xend, N_Vector y) {
  // reset counters
  nfcn = 0;
  njac = 0;
  nstep = 0;
  naccept = 0;
  nreject = 0;
  ndec = 0;
  nsol = 0;

  // Assuming autonomous (ifcn=0), explicit (imas=0), full Jacobian
  assert(deriv_fcn != NULL);
  assert(jac_fcn != NULL);
  assert(jac_decomp_fcn != NULL);
  assert(jac_solve_fcn != NULL);

  int k, kc, kopt;
  double hopt;

  create_work_vectors(y);
  std::vector<double> hh(km + 1, 0);
  std::vector<double> w(km + 1, 0);
  std::vector<double> a(km + 1, 0);
  std::vector<int> nj(km + 1, 0);

  setup_step_sequence(nj, a);
  x = *x_in;
  double xstart = x;

  double hmax_save = hmax;
  if (hmax == 0.0) hmax = xend - xstart;
  double posneg = xend - xstart > 0 ? 1.0 : -1.0;
  double min_rtol = N_VMin(rtol);
  double min_atol = N_VMin(atol);
  k = std::max(1,
               std::min(km - 1, int(-log10(min_rtol + min_atol) * 0.6 + 0.5)));
  hmaxn = std::min(fabs(hmax), fabs(xend - x));
  h = std::max(fabs(h), 1.0e-6);
  h = posneg * std::min(h, hmaxn);
  theta = 2 * fabs(thet);
  err = 0.0;

  /**/
  deriv_fcn(x, y, dy, user_data);
  if (output_fcn != NULL) {
    output_fcn(nstep, x, h, y, dy, output_fcn_data);
  }
  /**/

  w[0] = 1.0e30;
  N_VAbs(y, scal);
  N_VProd(rtol, scal, scal);
  N_VLinearSum(1.0, atol, 1.0, scal, scal);
  N_VInv(scal, scal);  // N.B. we work with inverse scal.

  caljac = false;
  reject = false;
  last = false;

g10:
  if (reject) theta = 2 * fabs(thet);
  atov = false;
  if (0.1 * fabs(xend - x) <= fabs(x) * uround) goto g110;
  hopt = h;
  h = posneg * std::min(std::min(fabs(h), fabs(xend - x)), hmaxn);
  if ((x + 1.01 * h - xend) * posneg > 0.0) {
    h = xend - x;
    last = true;
  }
  deriv_fcn(x, y, dy, user_data);
  nfcn += 1;

  // COMPUTATION OF THE JACOBIAN
  if (theta > thet && !caljac) {
    njac += 1;
    deriv_fcn(x, y, dy, user_data);
    jac_fcn(x, y, dy, user_data);
    caljac = true;
  }

  // FIRST AND LAST STEP
  if (nstep == 0 || last) {
    nstep += 1;
    for (int j = 0; j <= k; ++j) {
      kc = j;
      simex(x, y, j, hh, w, nj, a);
      if (atov) goto g10;
      if (j > 0 && err <= 1.0) goto g60;
    }
    goto g55;
  }

// --- BASIC INTEGRATION STEP
g30:
  nstep += 1;
  if (nstep >= nmax) goto g120;
  kc = k - 1;
  for (int j = 0; j <= kc; ++j) {
    simex(x, y, j, hh, w, nj, a);
    if (atov) goto g10;
  }

  // C *** *** *** *** *** *** ***
  // C --- CONVERGENCE MONITOR
  // C *** *** *** *** *** *** ***
  if (k == 1 || reject) goto g50;
  if (err <= 1.0) goto g60;
  if (err > pow(((double)(nj[k + 1] * nj[k])) / 4., 2.0)) goto g100;

g50:
  simex(x, y, k, hh, w, nj, a);
  if (atov) goto g10;
  kc = k;
  if (err <= 1.0) goto g60;

g55:
  if (err > pow(((double)nj[k + 1]) / 2., 2.0)) goto g100;
  kc = k + 1;
  // --- HOPE FOR CONVERGENCE IN LINE K+1
  simex(x, y, kc, hh, w, nj, a);
  if (atov) goto g10;
  if (err > 1.0) goto g100;

// C*** *** *** *** *** *** ***
// C --- STEP IS ACCEPTED
// C *** *** *** *** *** *** ***
g60:
  x = x + h;
  // TODO: Inefficient, also above
  N_VAbs(t[0], scal);  // make sure t[0] is set somewhere?
  N_VProd(rtol, scal, scal);
  N_VLinearSum(1.0, atol, 1.0, scal, scal);
  N_VInv(scal, scal);
  N_VScale(1.0, t[0], y);
  naccept += 1;
  caljac = false;
  // TODO: Interpolation functions for continous output
  if (output_fcn != NULL) {
    output_fcn(naccept, x, h, y, dy, output_fcn_data);
  }

  // C --- COMPUTE OPTIMAL ORDER
  if (kc == 1) {
    kopt = 2;
    if (reject) kopt = 1;
    goto g80;
  }
  if (kc <= k) {
    kopt = kc;
    if (w[kc - 1] < w[kc] * fac3) kopt = kc - 1;
    if (w[kc] < w[kc - 1] * fac4) kopt = std::min(kc + 1, km - 1);
  } else {
    kopt = kc - 1;
    if (kc > 2 && w[kc - 2] < w[kc - 1] * fac3) kopt = kc - 2;
    if (w[kc] < w[kopt] * fac4) kopt = std::min(kc, km - 1);
  }

g80:
  if (reject) {
    k = std::min(kopt, kc);
    h = posneg * std::min(fabs(h), fabs(hh[k]));
    reject = false;
    goto g10;
  }
  // C --- COMPUTE STEP SIZE FOR NEXT STEP
  if (kopt <= kc) {
    h = hh[kopt];
  } else {
    if (kc < k && w[kc] < w[kc - 1] * fac4) {
      h = hh[kc] * a[kopt + 1] / a[kc];
    } else {
      h = hh[kc] * a[kopt] / a[kc];
    }
  }
  k = kopt;
  h = posneg * fabs(h);
  goto g10;

// C *** *** *** *** *** *** ***
// C --- STEP IS REJECTED
// C *** *** *** *** *** *** ***
g100:
  k = std::min(k, kc);
  if (k > 1 && w[k - 1] < w[k] * fac3) k = k - 1;
  nreject += 1;
  h = posneg * hh[k];
  last = false;
  reject = true;
  if (caljac) goto g30;
  goto g10;

// C --- SOLUTION EXIT
g110:
  h = hopt;
  destroy_work_vectors();
  hmax = hmax_save;
  *x_in = x;
  return 0;

// C --- FAIL EXIT
g120:
  printf("  EXIT OF SODEX AT X=%g   H=%g  nstep=%d\n", x, h, nstep);
  destroy_work_vectors();
  hmax = hmax_save;
  *x_in = x;
  return 1;
}

// C *** *** *** *** *** *** ***
// C     S U B R O U T I N E    M I D E X
// C *** *** *** *** *** *** ***
int sodex::simex(double x, N_Vector y, int jj, std::vector<double>& hh,
                 std::vector<double>& w, std::vector<int>& nj,
                 std::vector<double>& a) {
  int m, njmid, ier;
  double fac, facmin, expo;
  // C --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE
  // C --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATE
  // C --- OF THE OPTIMAL STEP SIZE
  double hj = h / nj[jj];
  double hji = 1.0 / hj;

  // C *** *** *** *** *** *** ***
  // C  COMPUTE THE MATRIX E AND ITS DECOMPOSITION
  // C *** *** *** *** *** *** ***
  ier = jac_decomp_fcn(hj, user_data);
  ndec += 1;
  if (ier != 0) {
    goto g79;
  }
  // C *** *** *** *** *** *** ***
  // C --- STARTING PROCEDURE
  // C *** *** *** *** *** *** ***
  // ier = deriv_fcn(x+hj,y,dy,user_data);//?
  nfcn += 1;
  N_VScale(1.0, y, yh);
  N_VScale(1.0, dy, del);
  jac_solve_fcn(x, y, dy, del, tmp1, user_data);
  N_VScale(1.0, tmp1, del);
  N_VScale(hj, del, del);

  nsol += 1;
  m = nj[jj];

  // C *** *** *** *** *** *** ***
  // C --- SEMI-IMPLICIT MID-POINT METHOD
  // C *** *** *** *** *** *** ***
  // g30:
  for (int mm = 1; mm <= m; ++mm)  // was mm < m
  {
    N_VLinearSum(1.0, yh, 1.0, del, yh);

    ier = deriv_fcn(x + hj * mm, yh, dyh, user_data);
    if (ier != 0) goto g79;
    nfcn += 1;
    if (mm == m) continue;  // goto g30; //not sure about this?
    N_VLinearSum(2.0, dyh, -2.0 * hji, del, dyh);
    jac_solve_fcn(x, y, dyh, dyh, tmp1, user_data);
    nsol += 1;
    N_VScale(1.0, tmp1, dyh);
    N_VScale(hj, dyh, dyh);

    if (mm == 1 && jj <= 1)  // mm 1-based; jj 0-based
    {
      // C --- STABILITY CHECK
      double del1 = 0.0;
      del1 = N_VWL2Norm(del, scal);
      double del2 = 0.0;
      del2 = 0.5 * N_VWL2Norm(dyh, scal);
      theta = del2 / std::max(1.0, del1);
      if (theta > 1.0) goto g79;
    }
    N_VLinearSum(1.0, del, 1.0, dyh, del);
  }
  // C --- FINAL STEP (DUE TO BADER)
  N_VLinearSum(1.0, dyh, -hji, del, dyh);

  jac_solve_fcn(x, y, dyh, dyh, tmp1, user_data);
  N_VScale(1.0, tmp1, dyh);
  N_VScale(hj, dyh, dyh);
  nsol += 1;
  N_VLinearSum(1.0, yh, 1.0, dyh, t[jj]);  // t[jj] = yh + dyh

  // DEBUG N_VPrint_Serial(t[jj]);
  // C *** *** *** *** *** *** ***
  // C --- POLYNOMIAL EXTRAPOLATION
  // C *** *** *** *** *** *** ***
  if (jj == 0) return 0;
  for (int l = jj; l > 0; --l) {
    // double fac=((double) nj[jj])/((double)nj[l-1])-1.0; //TODO: OF tabulates
    // these
    double fac = pow(((double)nj[jj]) / ((double)nj[l - 1]), 2.0) - 1.0;
    fac = 1.0 / fac;
    N_VLinearSum(1.0, t[l], -1.0, t[l - 1], t[l - 1]);
    N_VLinearSum(1.0, t[l], fac, t[l - 1], t[l - 1]);
  }
  err = 0.0;
  // C --- SCALING
  N_VLinearSum(1.0, t[0], -1.0, t[1], wh);  // USING wh for temp storage here
  N_VAbs(wh, wh);
  err = N_VWrmsNorm(scal, wh);
  if (isnan(err)) goto g79;
  if (err >= 1.0e20) goto g79;
  if (err * uround >= 1.0) goto g79;
  if (jj > 1 && err >= errold) goto g79;
  errold = std::max(4 * err, 1.0);
  // C --- COMPUTE OPTIMAL STEP SIZES
  expo = 1.0 / (2.0 * (jj + 1) - 2);  // was 1/(2*jj-1) for odex
  facmin = pow(fac1, expo);
  fac = std::min(fac2 / facmin,
                 std::max(facmin, pow((err / safe1), expo) / safe2));
  fac = 1.0 / fac;
  hh[jj] = std::min(fabs(h) * fac, hmax);
  w[jj] = a[jj] / hh[jj];
  return 0;

g79:
  atov = true;
  h = h * 0.5;
  reject = true;
  return 1;
};

void sodex::setup_step_sequence(std::vector<int>& nj, std::vector<double>& a) {
  // --- DEFINE THE STEP SIZE SEQUENCE
  if (nsequ == 1) {
    nj[0] = 2;
    nj[1] = 6;
    if (km >= 2) nj[2] = 10;
    if (km >= 3) nj[3] = 14;
    if (km >= 4) nj[4] = 22;
    if (km >= 5) nj[5] = 34;
    if (km >= 6) nj[6] = 50;
    if (km >= 7) {
      for (int i = 7; i <= km; i++) {
        nj[i] = 2 * nj[i - 2] + nj[0];
      }
    }
  }

  a[0] = wkjac + (nj[0] + 1) * wkrow + wkdec;
  for (int i = 1; i <= km; i++) {
    a[i] = a[i - 1] + nj[i] * wkrow + wkdec;
  }
}

void sodex::create_work_vectors(N_Vector y) {
  yh = N_VClone(y);
  dy = N_VClone(y);
  dyh = N_VClone(y);
  del = N_VClone(y);
  wh = N_VClone(y);
  scal = N_VClone(y);
  tmp1 = N_VClone(y);
  tmp2 = N_VClone(y);
  tmp3 = N_VClone(y);
  t.resize(km + 1);
  for (int i = 0; i <= km; ++i) {
    t[i] = N_VClone(y);
  }
}

void sodex::destroy_work_vectors() {
  N_VDestroy(yh);
  N_VDestroy(dy);
  N_VDestroy(dyh);
  N_VDestroy(del);
  N_VDestroy(wh);
  N_VDestroy(scal);
  N_VDestroy(tmp1);
  N_VDestroy(tmp2);
  N_VDestroy(tmp3);
  for (int i = 0; i <= km; ++i) {
    N_VDestroy(t[i]);
  }
}

}  // end namespace sodex_cpp
