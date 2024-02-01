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
#ifndef OST_RANGE_HH
#define OST_RANGE_HH

#include <ost/module_config.hh> // for size_t, npos

namespace ost {
  
struct Range
{
  size_t Loc;
  size_t Length;
  Range()
    : Loc(String::npos), Length(0)
  { }

  Range(size_t loc,size_t len)
    : Loc(loc), Length(len)
  { }

  size_t End() const
  { return Loc+Length; }
  
  bool operator==(const Range& rhs) const 
  {
    return Loc==rhs.Loc && Length==rhs.Length;
  }
  bool operator!=(const Range& rhs) const 
  {
    return !this->operator==(rhs);
  }  
};

} // ost
#endif
