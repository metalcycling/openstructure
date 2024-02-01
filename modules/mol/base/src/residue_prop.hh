//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
#ifndef OST_RESIDUE_PROP_HH
#define OST_RESIDUE_PROP_HH

#include <vector>
#include <boost/operators.hpp>

#include <ost/mol/module_config.hh>

namespace ost { namespace mol {


class DLLEXPORT ResNum: private
    boost::additive<ResNum, int,
    boost::additive<ResNum,
    boost::totally_ordered<ResNum, 
    boost::totally_ordered<ResNum, int,
    boost::unit_steppable<ResNum> > > > >
{
public:
 
  // needed to wrap certain map classes
  ResNum():
    num_(1),alt_('\0')	 
  {}  

  ResNum(int n):
    num_(n), alt_('\0')
  { }

  ResNum(int n, char a):
    num_(n), alt_(a)
  {}

  bool operator==(const ResNum& r) const
  {
    return num_==r.num_ && alt_==r.alt_;
  }

  bool operator<(const ResNum& r) const
  {
    return num_==r.num_ ? alt_<r.alt_ : num_<r.num_;
  }

  int operator+=(int i)
  {
    num_+=i;
    return num_;
  }

  int operator-=(int i)
  {
    num_-=i;
    return num_;
  }

  int operator+=(const ResNum& r)
  {
    num_+=r.num_;
    return num_;
  }

  int operator-=(const ResNum& r)
  {
    num_-=r.num_;
    return num_;
  }
  
  ResNum& operator++()
  {
    ++num_;
    return *this;
  }

  ResNum& operator--()
  {
    --num_;
    return *this;
  }

  ResNum NextInsertionCode() const
  {
    char alt= alt_=='\0' ? 'a' : alt_+1;
    ResNum nrvo(num_,alt);
    return nrvo;
  }
  
  /// \brief get residue number as String
  ///
  /// The returned String consists of both the numeric residue sequence number
  /// and the insertion code. If the insertion code is not defined, i.e. is 
  /// equal to the null character, only the residue sequence number is returned.
  inline String AsString() const;
  
  int GetNum() const { return num_; }
  
  void SetNum(int num) { num_=num; }
  
  void SetInsCode(char ins_code) { alt_=ins_code; }
  
  char GetInsCode() const { return alt_; }  
  
private:
  int num_ : 24;
  int alt_ :  8;
};

typedef String ResidueKey;
typedef std::vector<ResNum> ResNumList;

inline std::ostream& operator<<(std::ostream& os, const ResNum& n)
{
  os << n.GetNum();
  if (n.GetInsCode()!='\0')
    os << n.GetInsCode();
  return os;
}

inline String ResNum::AsString() const
{
  std::stringstream ss;
  ss << *this;
  return ss.str();    
}

}} // ns


#endif
