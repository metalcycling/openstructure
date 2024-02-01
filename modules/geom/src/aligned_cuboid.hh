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
#ifndef GEOM_ALIGNED_CUBOID_HH
#define GEOM_ALIGNED_CUBOID_HH

/*
  Author: Marco Biasini
 */

#include "vec3.hh"
#include "module_config.hh"

namespace geom {


/// \brief axis-aligned cuboid
/// 
/// For an arbitrarily oriented cuboid see \ref Cuboid
class DLLEXPORT_OST_GEOM AlignedCuboid {
public:
  AlignedCuboid() :min_(), max_() {}
  AlignedCuboid(const Vec3& mmin, const Vec3& mmax) :min_(mmin), max_(mmax) {}
  
  Vec3 GetSize() const {return max_-min_;}
  
  Real GetVolume() const {
    Vec3 s=max_-min_;
    return s[0]*s[1]*s[2];
  }

  const Vec3& GetMin() const {return min_;}
  
  const Vec3& GetMax() const {return max_;}
  
  Vec3 GetCenter() const {return 0.5*(max_+min_);}
private:
  Vec3 min_;
  Vec3 max_;
};

inline std::ostream& operator<<(std::ostream& os, const AlignedCuboid& c)
{
  os << "(" << c.GetMin() << "," << c.GetMax() << ")";
  return os;
}

AlignedCuboid DLLEXPORT_OST_GEOM Union(const AlignedCuboid& lhs, const AlignedCuboid& rhs);

}


#endif
