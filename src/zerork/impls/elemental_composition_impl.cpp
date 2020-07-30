
#include <set>
#include <sstream> //ostringstream

#include "elemental_composition_impl.h"



namespace zerork {

ElementalComposition::Impl::Impl()
{
}

ElementalComposition::Impl::Impl(std::string name)
:
  name_(name)
{
}

ElementalComposition::Impl::Impl(const ElementalComposition::Impl& rhs)
:
  name_(rhs.name_),
  cmap_(rhs.cmap_)
{
}

ElementalComposition::Impl& ElementalComposition::Impl::operator=(Impl rhs)
{
  std::swap(this->name_,rhs.name_);
  std::swap(this->cmap_,rhs.cmap_);
  return *this;
}

int ElementalComposition::Impl::GetNumAtoms(const std::string& element) const
{
  int num_atoms = 0;
  CompositionMap::const_iterator iter = this->cmap_.find(element);
  if( iter != this->cmap_.end())
  {
      num_atoms = iter->second;
  }
  return num_atoms;
}


int ElementalComposition::Impl::GetNumHeavyAtoms() const
{
  int num_heavy_atoms = 0;
  for (CompositionMap::const_iterator it=this->cmap_.begin();
       it!=this->cmap_.end(); ++it)
  {
    if(it->first != std::string("H") && it->first != std::string("h"))
    {
      num_heavy_atoms += it->second;
    }
  }
  return num_heavy_atoms;
}

int ElementalComposition::Impl::GetNumTotalAtoms() const
{
  int num_total_atoms = 0;
  for (CompositionMap::const_iterator it=this->cmap_.begin();
       it!=this->cmap_.end(); ++it)
  {
    num_total_atoms += it->second;
  }
  return num_total_atoms;
}

std::vector<std::string> ElementalComposition::Impl::GetElementVector() const
{
  std::vector<std::string> elements;
  for (CompositionMap::const_iterator it=this->cmap_.begin();
       it!=this->cmap_.end(); ++it)
  {
    if(it->second > 0) elements.push_back(it->first);
  }
  return elements;
}

std::string& ElementalComposition::Impl::name()
{
  return this->name_;
}

const std::string& ElementalComposition::Impl::name() const
{
  return this->name_;
}

std::string ElementalComposition::Impl::ToString() const
{
  std::ostringstream oss;
  for (CompositionMap::const_iterator it=this->cmap_.begin();
       it!=this->cmap_.end(); ++it)
  {
    if(it->second > 0) {
      oss << it->first;
    }
    if(it->second > 1) {
      oss << it->second;
    }
  }
  return oss.str();
}

std::string ElementalComposition::Impl::ToStringWithSeparator(std::string sep) const
{
  std::ostringstream oss;
  for (CompositionMap::const_iterator it=this->cmap_.begin();
       it!=this->cmap_.end(); ++it)
  {
    if(it->second > 0) {
      oss << it->first << sep << it->second << " ";
    }
  }
  std::string oss_str = oss.str();
  return oss_str.substr(0,oss_str.size()-1); //delete last space
}

int& ElementalComposition::Impl::operator[](const std::string& element)
{
  return this->cmap_[element];
}

int ElementalComposition::Impl::operator[](const std::string& element) const
{
  CompositionMap::const_iterator iter = this->cmap_.find(element);
  if(iter != this->cmap_.end())
  {
      return iter->second;
  }
  return 0;
}

void ElementalComposition::Impl::clear()
{
  this->cmap_.clear();
}

int ElementalComposition::Impl::compare(const Impl& rhs) const
{
  //{ /* do actual comparison */ }
  const Impl& lhs = *this;

  int lsize = lhs.GetNumTotalAtoms();
  int rsize = rhs.GetNumTotalAtoms();
  if(lsize > rsize) {
    return +1; // lhs goes second if more elements
  } else if(lsize < rsize) {
    return -1;  // rhs goes second if more elements
  } else {
    std::set<std::string> elements;
    //Iterate over elements of l
    for (CompositionMap::const_iterator it=lhs.cmap_.begin(); it!=lhs.cmap_.end(); ++it)
    {
      elements.insert(it->first);
    }
    //Iterate over elements of r
    for (CompositionMap::const_iterator it=rhs.cmap_.begin(); it!=rhs.cmap_.end(); ++it)
    {
      elements.insert(it->first);
    }

    for (std::set<std::string>::const_iterator it=elements.begin(); it!=elements.end(); ++it)
    {
      std::string e = *it;
      int lval = 0;
      int rval = 0;
      CompositionMap::const_iterator l_iter = lhs.cmap_.find(e);
      if( l_iter != lhs.cmap_.end())
      {
         lval = l_iter->second;
      }
      CompositionMap::const_iterator r_iter = rhs.cmap_.find(e);
      if( r_iter != rhs.cmap_.end())
      {
         rval = r_iter->second;
      }
      if(lval > rval) return +1;
      if(lval < rval) return -1;
    }
    if(lhs.name_ > rhs.name_) return +1;
    if(lhs.name_ < rhs.name_) return -1;
  }
  return 0;
}


bool ElementalComposition::Impl::operator==(const ElementalComposition::Impl& rhs) const
{
  return this->compare(rhs)==0;
}

bool ElementalComposition::Impl::operator!=(const ElementalComposition::Impl& rhs) const
{
  return !this->operator==(rhs);
}

bool ElementalComposition::Impl::operator< (const ElementalComposition::Impl& rhs) const
{
  return this->compare(rhs)==-1;
}

bool ElementalComposition::Impl::operator> (const ElementalComposition::Impl& rhs) const
{
  return  rhs.operator< (*this);
}

bool ElementalComposition::Impl::operator<=(const ElementalComposition::Impl& rhs) const
{
  return !this->operator> (rhs);
}

bool ElementalComposition::Impl::operator>=(const ElementalComposition::Impl& rhs) const
{
  return !this->operator< (rhs);
}


//Arithmetic operators
ElementalComposition::Impl ElementalComposition::Impl::operator+(const ElementalComposition::Impl& rhs)
{
  Impl tmp(*this);
  tmp += rhs;
  return tmp;
}

ElementalComposition::Impl ElementalComposition::Impl::operator-(const ElementalComposition::Impl& rhs)
{
  Impl tmp(*this);
  tmp -= rhs;
  return tmp;
}

ElementalComposition::Impl& ElementalComposition::Impl::operator+=(const ElementalComposition::Impl& rhs)
{
  for(CompositionMap::const_iterator it=rhs.cmap_.begin();
       it!=rhs.cmap_.end(); ++it)
  {
    this->cmap_[it->first] += it->second;
  }
  return *this;
}

ElementalComposition::Impl& ElementalComposition::Impl::operator-=(const ElementalComposition::Impl& rhs)
{
  for(CompositionMap::const_iterator it=rhs.cmap_.begin();
       it!=rhs.cmap_.end(); ++it)
  {
    this->cmap_[it->first] -= it->second;
  }
  return *this;
}

} //end namespace zerork
