
#include "elemental_composition.h"
#include "impls/elemental_composition_impl.h"


namespace zerork {


ElementalComposition::ElementalComposition()
:
  impl_(new Impl())
{
}

ElementalComposition::ElementalComposition(std::string name)
:
  impl_(new Impl(name))
{
}


ElementalComposition::ElementalComposition(const ElementalComposition& rhs)
:
  impl_(new Impl(*rhs.impl_))
{
}

ElementalComposition::~ElementalComposition()
{
  //unique_ptr deletes impl_
}

ElementalComposition& ElementalComposition::operator=(ElementalComposition rhs)
{
  std::swap(this->impl_, rhs.impl_);
  return *this;
}

//Getters
int ElementalComposition::GetNumAtoms(const std::string& element) const
{
  return impl_->GetNumAtoms(element);
}
int ElementalComposition::GetNumHeavyAtoms() const
{
  return impl_->GetNumHeavyAtoms();
}
int ElementalComposition::GetNumTotalAtoms() const
{
  return impl_->GetNumTotalAtoms();
}

std::vector<std::string> ElementalComposition::GetElementVector() const
{
  return impl_->GetElementVector();
}

std::string& ElementalComposition::name()
{
  return impl_->name();
}

const std::string& ElementalComposition::name() const
{
  return impl_->name();
}

std::string ElementalComposition::ToString() const
{
  return impl_->ToString();
}

std::string ElementalComposition::ToStringWithSeparator(std::string sep) const
{
  return impl_->ToStringWithSeparator(sep);
}

int& ElementalComposition::operator[](const std::string& element)
{
  return impl_->operator[](element);
}

const int& ElementalComposition::operator[](const std::string& element) const
{
  return impl_->operator[](element);
}

void ElementalComposition::clear()
{
  impl_->clear();
}

bool operator==(const ElementalComposition& lhs, const ElementalComposition& rhs)
{
  return lhs.impl_->operator==(*rhs.impl_);
}

bool operator!=(const ElementalComposition& lhs, const ElementalComposition& rhs)
{
  return lhs.impl_->operator!=(*rhs.impl_);
}

bool operator< (const ElementalComposition& lhs, const ElementalComposition& rhs)
{
  return lhs.impl_->operator< (*rhs.impl_);
}

bool operator> (const ElementalComposition& lhs, const ElementalComposition& rhs)
{
  return lhs.impl_->operator> (*rhs.impl_);
}

bool operator<=(const ElementalComposition& lhs, const ElementalComposition& rhs)
{
  return lhs.impl_->operator<=(*rhs.impl_);
}

bool operator>=(const ElementalComposition& lhs, const ElementalComposition& rhs)
{
  return lhs.impl_->operator>=(*rhs.impl_);
}

ElementalComposition operator+(ElementalComposition lhs, const ElementalComposition& rhs)
{
  *(lhs.impl_) += *(rhs.impl_);
  return lhs;
}

ElementalComposition operator-(ElementalComposition lhs, const ElementalComposition& rhs)
{
  *(lhs.impl_) -= *(rhs.impl_);
  return lhs;
}

ElementalComposition& ElementalComposition::operator+=(const ElementalComposition& rhs)
{
  *(this->impl_) += *(rhs.impl_);
  return *this;
}

ElementalComposition& ElementalComposition::operator-=(const ElementalComposition& rhs)
{
  *(this->impl_) -= *(rhs.impl_);
  return *this;
}

} //end namespace zerork
