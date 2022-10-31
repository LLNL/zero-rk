#ifndef NASA_POLY_H
#define NASA_POLY_H

namespace zerork {

class nasa_poly_group
{
 public:
  nasa_poly_group(const int inpSpc, double const * const inpCoef,
		  double const * const inpTlow, double const * const inpThigh);

  virtual ~nasa_poly_group();
  void getCp_R(const double T, double Cp_R[]) const;
  void getH_RT(const double T, double H_RT[]) const;
  void getG_RT(const double T, double G_RT[]) const;

  void getCp_R_mr(const int nReactors, const double T[], double Cp_R[]) const;
  void getH_RT_mr(const int nReactors, const double T[], double H_RT[]) const;

  void getThermoCoeffs(double coeffs[]) const;

 protected:
  int nGroupSpc;
  double *thermoCoef;
  double *Tlow;
  double *Thigh;
  // the enthalpy coefficient multipliers are padded with ones at the
  // beginning as in the gpu implementation
};

} // namespace zerork


#endif
