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
#ifndef OST_MOL_ALG_FIND_MEMBRANE_H
#define OST_MOL_ALG_FIND_MEMBRANE_H

#include <ost/mol/mol.hh>
#include <ost/geom/geom.hh>

namespace ost { namespace mol { namespace alg {

struct FindMemParam{
  FindMemParam() { }

  geom::Vec3 GetMembraneAxis() const;
  geom::Vec3 axis;
  geom::Vec3 tilt_axis;
  Real tilt;
  Real angle;
  Real width;
  Real pos;
  Real energy;
  Real membrane_asa;
  ost::mol::EntityHandle membrane_representation;
};

FindMemParam FindMembrane(ost::mol::EntityHandle& ent, 
                          bool assign_membrane_representation,
                          bool fast);

FindMemParam FindMembrane(ost::mol::EntityView& ent, 
                          bool assign_membrane_representation,
                          bool fast);

}}} // ns

#endif
