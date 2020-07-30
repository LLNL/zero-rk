#ifndef REACTOR_CONSTANTS_H_
#define REACTOR_CONSTANTS_H_

#include "stdio.h" // needed for printf

enum MatrixType {DENSE_COL_MAJOR, COMPRESSED_COL_STORAGE};

enum ReactorError {NONE,
                   INDEX_OUT_OF_RANGE,
                   INVALID_MECHANISM,
                   UNKNOWN};

// ---------------------------------------------------------------------------
// Abstract base class for the heat loss function. Users are able to define
// a derived HeatLoss class by defining the constructor, destrutor and
// GetHeatLoss function according to the pure virtual prototype below.
class HeatLoss
{
 public:
  virtual ~HeatLoss() {};
  virtual ReactorError GetHeatLoss(const double reactor_time,
                                   const double state[],
                                   double *heat_loss) = 0;
  // user is responsible for defining the GetHeatLoss function that returns 
  // the heat_loss [J/s] (via pointer) computed at a given reactor_time,
  // and given system state.  The function returns a ReactorError 
  // value to signal to the caller if an error occured.
};

// Derived HeatLoss class for an adiabatic (zero) heat loss reactor
class Adiabatic : public HeatLoss
{
 public:
  Adiabatic() {};   // no parameters or memory needs to be set
  ~Adiabatic() {};  // no parameters or memory needs to be freed
  ReactorError GetHeatLoss(const double reactor_time,
                           const double state[],
                           double *heat_loss) {
    *heat_loss = 0.0;
    return NONE;      // no reactor error
  }
};
// ---------------------------------------------------------------------------
// Abstract base class for the volume function. Users are able to define
// a derived Volume class by defining the constructor, destrutor and
// GetVolume function according to the pure virtual prototype below.
class Volume
{
 public:
  virtual ~Volume() {};
  virtual ReactorError GetVolume(const double reactor_time,
                                 const double state[],
                                 double *volume,
                                 double *dvolume_dt) = 0;
  // user is responsible for defining the GetHeatLoss function that returns 
  // the volume [m^3] (via pointer), and time derivative of volume [m^3/s]
  // (via dvolume_dt pointer) computed at a given reactor_time, and given
  // system state.  The function returns a ReactorError value 
  // to signal to the caller if an error occured.
};

// Derived volume class for a constant volume reactor with unit volume 1 [m^3]
class ConstantUnitVolume : public Volume
{
 public:
  ConstantUnitVolume() {};  // no parameters or memory needs to be set
  ~ConstantUnitVolume() {}; // no parameters or memory needs to be freed
  ReactorError GetVolume(const double reactor_time,
                         const double state[],
                         double *volume,
                         double *dvolume_dt) {
    *volume     = 1.0;
    *dvolume_dt = 0.0;
    return NONE;      // no reactor error
  }
};

// ---------------------------------------------------------------------------
                     
// header functions must be inlined
inline bool CheckReactorError(const char context[],
                              const ReactorError error)
{
  char preface[] = "ReturnedError:";

  switch(error) {

  case NONE:
    return false;

  case INDEX_OUT_OF_RANGE:
    printf("%s %s [%s]\n",preface,"index out of range",context);
    return true;

  case INVALID_MECHANISM:
    printf("%s %s [%s]\n",preface,"invalid mechanism",context);
    return true;

  case UNKNOWN:
    printf("%s %s [%s]\n",preface,"unknown error",context);
    return true;

  default:
    printf("%s %s [%s]\n",preface,"undefined error",context);
    return true;
  }
}

#endif
