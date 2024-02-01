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
#include "entity_property_mapper.hh"
#include <ost/generic_property.hh>
/*
  Author: Marco Biasini
 */

namespace ost { namespace mol {

EntityPropertyMapper::EntityPropertyMapper(const String& prop_name, 
                                           Prop::Level level):
  prop_(PropertyFromString(prop_name)), prop_name_(prop_name)
{
  if (prop_.id==Prop::UNDEF) {
    prop_.level=level;
  }
}

template <typename T, bool B>
Real EntityPropertyMapper::get_property(const T& atom, Real def_value) const
{
  if (prop_.id!=Prop::UNDEF) {
    switch (prop_.level) {
      case Prop::ATOM:
          return atom.GetFloatProperty(prop_.id);
      default:
        switch(prop_.type) {
          case Prop::FLOAT:
            return atom.GetResidue().GetFloatProperty(prop_.id);
          default:
            return atom.GetResidue().GetIntProperty(prop_.id);
        }
    }
  }
  switch(prop_.level) {
    case Prop::ATOM:
      if (B) {
        return atom.GetFloatProp(prop_name_, def_value);
      } else {
        return atom.GetFloatProp(prop_name_);          
      }
    case Prop::RESIDUE:
      if (B) {
        Real value=atom.GetResidue().GetFloatProp(prop_name_, 
                                                             def_value);
        return value; 
      } else {
        return atom.GetResidue().GetFloatProp(prop_name_);
      }
    default:
     return def_value;
  }

}

Real EntityPropertyMapper::Get(const AtomView& atom) const
{
  return this->get_property<AtomView, false>(atom, 0.0);
}

Real EntityPropertyMapper::Get(const AtomHandle& atom) const
{
  return this->get_property<AtomHandle, false>(atom, 0.0);
}



Real EntityPropertyMapper::Get(const AtomHandle& atom, 
                                 Real def_value) const
{
  return this->get_property<AtomHandle, true>(atom, def_value);
}

Real EntityPropertyMapper::Get(const AtomView& atom, 
                                 Real def_value) const
{
  return this->get_property<AtomView, true>(atom, def_value);
}



Real EntityPropertyMapper::Get(const ResidueView& res) const
{
  return res.GetFloatProp(prop_name_);
}

Real EntityPropertyMapper::Get(const ResidueHandle& res) const
{
  return res.GetFloatProp(prop_name_);
}



Real EntityPropertyMapper::Get(const ResidueHandle& res,
                                 Real def_value) const
{
  return res.GetFloatProp(prop_name_, def_value);
}

Real EntityPropertyMapper::Get(const ResidueView& res,
                                 Real def_value) const
{
  return res.GetFloatProp(prop_name_, def_value);
}






Real EntityPropertyMapper::Get(const ChainView& chain) const
{
  return chain.GetFloatProp(prop_name_);
}

Real EntityPropertyMapper::Get(const ChainHandle& chain) const
{
  return chain.GetFloatProp(prop_name_);
}



Real EntityPropertyMapper::Get(const ChainHandle& chain,
                                 Real def_value) const
{
  return chain.GetFloatProp(prop_name_, def_value);
}

Real EntityPropertyMapper::Get(const ChainView& chain,
                                 Real def_value) const
{
  return chain.GetFloatProp(prop_name_, def_value);
}

}} //ns
