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
#ifndef OST_MOL_BOUNDING_BOX_HH
#define OST_MOL_BOUNDING_BOX_HH

#include <ost/geom/geom.hh>
#include <ost/mol/module_config.hh>
/*
  Author: Marco Biasini
 */
namespace ost { namespace mol {

class EntityHandle;
class EntityView;

/// \name create bounding box from entity
//@{
/// \brief create bounding box from entity handle
/// \sa \ref bbox.py "Bounding Box Example"
geom::Cuboid DLLEXPORT_OST_MOL BoundingBoxFromEntity(const EntityHandle& ent);

/// \brief create bounding box from entity view
/// \sa \ref bbox.py "Bounding Box Example"
geom::Cuboid DLLEXPORT_OST_MOL BoundingBoxFromEntity(const EntityView& ent);

/// \brief create bounding box from atom list
/// \sa \ref bbox.py "Bounding Box Example"
geom::Cuboid DLLEXPORT_OST_MOL BoundingBoxFromAtoms(const AtomHandleList& atoms);
//@}

/// \example bbox.py
/// Illustrates the use of BoundingBoxFromEntity().
}}

#endif
