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
#ifndef OST_ATOM_PROP_HH
#define OST_ATOM_PROP_HH

#include <map>
#include <ost/mol/module_config.hh>
#include <ost/geom/mat3.hh>

namespace ost { namespace mol { namespace impl {


class DLLEXPORT_OST_MOL AtomProp {
public:
  AtomProp()
    : radius(0.0), charge(0.0), mass(0.0), has_anisou(false), is_default(true),
      anisou(geom::Mat3()) 
  {
  }
  AtomProp(Real r, Real m, Real c=0.0): radius(r), charge(c), 
     mass(m), has_anisou(false), is_default(true)
  { }
  AtomProp(const AtomProp& rhs): radius(rhs.radius), charge(rhs.charge),
    mass(rhs.mass), has_anisou(rhs.has_anisou), is_default(rhs.is_default),
    anisou(rhs.anisou)
  { }
  Real          radius;
  Real          charge;
  Real          mass;
  bool          has_anisou;
  bool          is_default;
  geom::Mat3    anisou;
  
  static AtomProp* GetDefaultProps(const String& ele);
  
private:
  static void Init();
  static std::map<String, AtomProp> defaults_;
};

}}} // ns

#endif

