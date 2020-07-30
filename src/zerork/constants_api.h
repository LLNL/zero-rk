#ifndef ZERORK_CONSTANTS_API_H_
#define ZERORK_CONSTANTS_API_H_

#include <stddef.h> // define size_t

namespace zerork {

class PhysicalConstants
{
 public:
  PhysicalConstants();
  ~PhysicalConstants();
  
  double GetGasConstant() const;                 // [J/kmol/K]
  double GetBoltzmannConstant() const;           // [J/K]
  double GetAvogadroNum() const;                 // [particles/kmol]
  double GetJoulesPerCalorie() const;            // [J/cal]
  double GetAtomicMass(size_t atomic_num) const; // [kg/kmol] If the atomic 
  // number is not in the range [1,92] inclusively, then zero is returned.
                                           
 private:
  PhysicalConstants(const PhysicalConstants&);            // noncopyable
  PhysicalConstants& operator=(const PhysicalConstants&); // noncopyable
  class Impl;
  Impl *impl_;
};

class ZeroRKConstants
{
 public:
  ZeroRKConstants();
  ~ZeroRKConstants();

  // The commit id is from the shell command `git rev-parse HEAD`
  const char * GetCommitId() const; 
  // The commit timestamp is from the shell command `git show -s --format=%ci`
  const char * GetCommitTimestamp() const;
  // The commit branch is from the shell command 
  //   `git rev-parse --abbrev-ref HEAD`
  const char * GetBranchName() const;
  // The type of fast exponential
  const char * GetExpType() const;

 private:
  ZeroRKConstants(const ZeroRKConstants&);            // noncopyable
  ZeroRKConstants& operator=(const ZeroRKConstants&); // noncopyable
  class Impl;
  Impl *impl_;
};


} //end namespace zerork

#endif
