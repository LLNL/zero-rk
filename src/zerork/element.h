#ifndef ZERORK_ELEMENT_H
#define ZERORK_ELEMENT_H

namespace zerork {

class element
{
 public:

  element();
  element(const int num, const double m, const char *sym);
  element(const int num);
 
  ~element();

  // data assignment
  bool setElement(const int num, const double m, const char *sym);
  bool setElement(const int num);
  bool setElement(const char *sym);

  // data access (constant member functions)
  int    getNumber() const {return number;}
  double getMass() const {return mass;}
  const char * getSymbol_c_str() const {return symbol;}
  void   printInfo() const;

 private:
  int    number;
  double mass;
  char   symbol[3];
};

} // namespace zerork

#endif
