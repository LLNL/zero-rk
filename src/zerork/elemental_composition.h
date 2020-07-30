#ifndef ELEMENTAL_COMPOSITION_H_
#define ELEMENTAL_COMPOSITION_H_

#include <string>
#include <vector>
#include <memory> //unique_ptr

namespace zerork {

class ElementalComposition
{
 public:
  ElementalComposition();
  ElementalComposition(std::string name);
  ElementalComposition(const ElementalComposition&);
  ElementalComposition& operator=(ElementalComposition rhs);
  virtual ~ElementalComposition();

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
  const int& operator[](const std::string&) const;
  void clear();

  //Comparison operators
  friend bool operator==(const ElementalComposition& lhs, const ElementalComposition& rhs);
  friend bool operator!=(const ElementalComposition& lhs, const ElementalComposition& rhs);
  friend bool operator< (const ElementalComposition& lhs, const ElementalComposition& rhs);
  friend bool operator> (const ElementalComposition& lhs, const ElementalComposition& rhs);
  friend bool operator<=(const ElementalComposition& lhs, const ElementalComposition& rhs);
  friend bool operator>=(const ElementalComposition& lhs, const ElementalComposition& rhs);

  //Arithmetic operators
  friend ElementalComposition operator+(ElementalComposition lhs, const ElementalComposition& rhs);
  friend ElementalComposition operator-(ElementalComposition lhs, const ElementalComposition& rhs);
  ElementalComposition& operator+=(const ElementalComposition& rhs);
  ElementalComposition& operator-=(const ElementalComposition& rhs);


 private:
  class Impl;
  std::unique_ptr<ElementalComposition::Impl> impl_;
};

//Comparison operators
bool operator==(const ElementalComposition& lhs, const ElementalComposition& rhs);
bool operator!=(const ElementalComposition& lhs, const ElementalComposition& rhs);
bool operator< (const ElementalComposition& lhs, const ElementalComposition& rhs);
bool operator> (const ElementalComposition& lhs, const ElementalComposition& rhs);
bool operator<=(const ElementalComposition& lhs, const ElementalComposition& rhs);
bool operator>=(const ElementalComposition& lhs, const ElementalComposition& rhs);

//Binary arithmetic operators
ElementalComposition operator+(ElementalComposition lhs, const ElementalComposition& rhs);
ElementalComposition operator-(ElementalComposition lhs, const ElementalComposition& rhs);

} //end namespace zerork

#endif
