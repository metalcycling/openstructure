//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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

/*
  convenience null function

  Author: Ansgar Philippsen
*/

#ifndef IPLT_NULL_FUNCTION_H
#define IPLT_NULL_FUNCTION_H

#include "function_base.hh"

namespace ost { namespace iplt {

//! null function
/*!
  Implements Function interface, will always
  return zero
*/
class DLLEXPORT_OST_IPLT_BASE NullFunction: public Function {
public:
  NullFunction():
    Function(SPATIAL) 
  {}
  
  virtual DataType GetType() const {return REAL;}

  virtual Real GetIntpolReal(const Vec3& v) const {return 0.0;}

  virtual Complex GetIntpolComplex(const Vec3& v) const {return Complex(0.0,0.0);}
};

}} // namespace


#endif
