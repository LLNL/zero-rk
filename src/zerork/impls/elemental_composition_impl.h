#ifndef ELEMENTAL_COMPOSITION_IMPL_H_
#define ELEMENTAL_COMPOSITION_IMPL_H_

#include "../elemental_composition.h"

#include <map>
namespace zerork {


class ElementalComposition::Impl
{
 public:
  Impl();
  Impl(std::string name);
  Impl(const Impl&); //Copy constructor
  Impl& operator=(Impl rhs);
  virtual ~Impl() {};

  //Getters
  int GetNumAtoms(const std::string&) const;
  int GetNumHeavyAtoms() const;
  int GetNumTotalAtoms() const;
  std::vector<std::string> GetElementVector() const;
  std::string& name();
  const std::string& name() const;
  std::string ToString() const;
  std::string ToStringWithSeparator(std::string) const;

  //Access operators
  int& operator[](const std::string&);
  int operator[](const std::string&) const;
  void clear();

  //Comparison operators
  //friend bool operator==(const Impl& lhs, const Impl& rhs);
  //friend bool operator!=(const Impl& lhs, const Impl& rhs);
  //friend bool operator< (const Impl& lhs, const Impl& rhs);
  //friend bool operator> (const Impl& lhs, const Impl& rhs);
  //friend bool operator<=(const Impl& lhs, const Impl& rhs);
  //friend bool operator>=(const Impl& lhs, const Impl& rhs);

  int compare(const Impl& rhs) const;

  bool operator==(const Impl& rhs) const;
  bool operator!=(const Impl& rhs) const;
  bool operator< (const Impl& rhs) const;
  bool operator> (const Impl& rhs) const;
  bool operator<=(const Impl& rhs) const;
  bool operator>=(const Impl& rhs) const;

  //Arithmetic operators
  //friend Impl operator+(Impl lhs, const Impl& rhs);
  //friend Impl operator-(Impl lhs, const Impl& rhs);
  Impl operator+(const Impl& rhs);
  Impl operator-(const Impl& rhs);
  Impl& operator+=(const Impl& rhs);
  Impl& operator-=(const Impl& rhs);

 private:
  std::string name_;
  typedef std::map<std::string,int> CompositionMap;
  // Consider custom sort on cmap to favor heavy atoms
  // in comparison operators
  CompositionMap cmap_;
};

//N.B. Can't define as binary functions because Impl is private.
//bool operator==(const ElementalComposition::Impl& lhs, const ElementalComposition::Impl& rhs);
//bool operator!=(const ElementalComposition::Impl& lhs, const ElementalComposition::Impl& rhs);
//bool operator< (const ElementalComposition::Impl& lhs, const ElementalComposition::Impl& rhs);
//bool operator> (const ElementalComposition::Impl& lhs, const ElementalComposition::Impl& rhs);
//bool operator<=(const ElementalComposition::Impl& lhs, const ElementalComposition::Impl& rhs);
//bool operator>=(const ElementalComposition::Impl& lhs, const ElementalComposition::Impl& rhs);

//ElementalComposition::Impl operator+(ElementalComposition::Impl lhs, const ElementalComposition::Impl& rhs);
//ElementalComposition::Impl operator-(ElementalComposition::Impl lhs, const ElementalComposition::Impl& rhs);

} // end namespace zerork

#endif
