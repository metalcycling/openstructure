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
#include "mat.hh"
#include "exc.hh"

#include "vecmat3_op.hh"
#include "vec3.hh"

namespace geom {

geom::Vec3 Mat3::GetCol(int index) const
{
  return geom::Vec3(data_[0][index], data_[1][index], data_[2][index]);  
}

geom::Vec3 Mat3::GetRow(int index) const
{
  return geom::Vec3(data_[index][0], data_[index][1], data_[index][2]);
}


////////////////////////////////////////////
// private member functions

std::ostream& operator<<(std::ostream& os, const Mat3& m)
{
  os << "{{" << m(0,0) << "," << m(0,1) << "," << m(0,2) << "},";
  os << "{" << m(1,0) << "," << m(1,1) << "," << m(1,2) << "},";
  os << "{" << m(2,0) << "," << m(2,1) << "," << m(2,2) << "}}";
  return os;
}

} // ns
