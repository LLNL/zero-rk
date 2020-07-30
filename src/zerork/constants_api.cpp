#include "constants.h"
#include "atomicMassDB.h"

#include "constants_api.h"

#include "zerork_version_info.h"
#include "fast_exps.h"

// Make the git version information available to this source file.
namespace zerork {

// ---------------------------------------------------------------------------
class PhysicalConstants::Impl
{
 public:
  Impl() {};
  double GetGasConstant() const {return zerork::NIST_RU;} 
  double GetBoltzmannConstant() const {return zerork::KBoltzmann;}
  double GetAvogadroNum() const {return zerork::KAvogadroNumber;}
  double GetJoulesPerCalorie() const 
    {return zerork::CAL_PER_MOL_TACT*zerork::NIST_RU;}
  double GetAtomicMass(size_t atomic_num) const;
};

double PhysicalConstants::Impl::GetAtomicMass(size_t atomic_num) const
{
  if(1 <= atomic_num && atomic_num <= zerork::MAX_ELEMENTS) {
    return zerork::atomicMassDB[atomic_num-1];
  }
  // 
  return 0.0;
}

// ---------------------------------------------------------------------------
class ZeroRKConstants::Impl
{
 public:
  Impl() {};

  const char * GetCommitId() const {return zerork::git_commit_id;}
  const char * GetCommitTimestamp() const {return zerork::git_commit_timestamp;}
  const char * GetBranchName() const {return zerork::git_branch_name;}
  const char * GetExpType() const {return zerork::expType;}
};

// ---------------------------------------------------------------------------
// Public facing API for the PhysicalConstants class
PhysicalConstants::PhysicalConstants()
{
  impl_ = new Impl();
}

PhysicalConstants::~PhysicalConstants()
{
  if(impl_ != NULL) {
    delete impl_;
  }
}

double PhysicalConstants::GetGasConstant() const
{
  return impl_->GetGasConstant();
}

double PhysicalConstants::GetBoltzmannConstant() const
{
  return impl_->GetBoltzmannConstant();
}

double PhysicalConstants::GetAvogadroNum() const
{
  return impl_->GetAvogadroNum();
}

double PhysicalConstants::GetJoulesPerCalorie() const
{
  return impl_->GetJoulesPerCalorie();
}

double PhysicalConstants::GetAtomicMass(size_t atomic_num) const
{
  return impl_->GetAtomicMass(atomic_num);
}

// ---------------------------------------------------------------------------
// Public facing API for the ZeroRKConstants class
ZeroRKConstants::ZeroRKConstants()
{
  impl_ = new Impl();
}

ZeroRKConstants::~ZeroRKConstants()
{
  if(impl_ != NULL) {
    delete impl_;
  }
}

const char * ZeroRKConstants::GetCommitId() const 
{
  return impl_->GetCommitId();
}

const char * ZeroRKConstants::GetCommitTimestamp() const 
{
  return impl_->GetCommitTimestamp();
}

const char * ZeroRKConstants::GetBranchName() const 
{
  return impl_->GetBranchName();
}

const char * ZeroRKConstants::GetExpType() const 
{
  return impl_->GetExpType();
}

} // namespace zerork

