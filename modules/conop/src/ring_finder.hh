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
#ifndef OST_CONOP_RING_FINDER_HH
#define OST_CONOP_RING_FINDER_HH

/*
   Author: Tobias Schmidt
 */

#include <ost/conop/module_config.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/entity_view.hh>
#include <map>

namespace ost { namespace conop {

/// \brief  Class for finding rings in an entity
///
/// This is an implementation of a depth first search based method to identify
/// all atoms and bonds belonging to rings in an entity. The time complexity is
/// linear. The search for rings is only carried out once and the ring
/// information is stored in generic properties in the entity and it can be
/// obtained as an EntityView.
class DLLEXPORT_OST_CONOP RingFinder {
public:
  RingFinder(mol::EntityHandle &eh);
public:

  /// \brief Walk along the path of the entity in a depth first search and
  ///        identify atoms and bonds which are in rings
  ///
  /// For all atoms and bonds, belonging to a ring, the BoolProp
  /// RFinring is set to true. The search for rings is only carried out once.
  /// This method does not need to be called explicitly, it will be
  /// automatically called if the rings have not yet been identified.
  void PerceiveRings();

  /// \brief Get an EntityView, containing all atoms and bonds belonging to
  ///        rings
  mol::EntityView GetRings();

  /// \brief Get number of atoms belonging to rings
  int GetRingAtomCount();

  /// \brief Get number of bonds belonging to rings
  int GetRingBondCount();

  /// \brief Check if rings have been perceived
  bool RingsPerceived();

private:
  int VisitNext(mol::AtomHandle &curr, mol::AtomHandle &prev,
                std::map<long int,bool> &ring_closings);
  mol::EntityHandle ent_;
  mol::EntityView ring_view_;
  bool rings_perceived_;
};

}}

#endif
