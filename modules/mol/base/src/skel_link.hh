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
#ifndef OST_SKEL_LINK_HH
#define OST_SKEL_LINK_HH

#include "skel_link_fw.hh"
#include "atom_fw.hh"

namespace ost { namespace mol {

/*
  Provides low-level skeleton link connectivity

  xi denotes rotation around a link, which will only affect the
  consecutive ones, ie the direction normal. In a molecular structure
  xi denotes a torsion angle.

  omega and xi denote orientation relativ to the current transformation,
  which has a main axis direction and a normal vector. If the main axis
  is along z, and the normal vector along x, omega describes the angle
  of the link with the z-axis, and xi the angle of the projected link
  in the xy plane to the x-axis.
*/

class DLLEXPORT_OST_MOL SkelLink {
public:
  // initialize with only two atoms
  /*
    since only the length of the link can be determined, the
    angles are set to zero, and a call to RebuildICS() is
    required to update these.
  */
  SkelLink(const AtomHandle& a1, const AtomHandle& a2);

  // initialize with two atoms and angle parameters
  /*
   */
  SkelLink(const AtomHandle& a1, const AtomHandle& a2, Real len, Real xi, Real omega, Real xi);

  // rebuild ics from xcs
  void RebuildICS();

  // rebuild xcs from ics
  void RebuildXCS();

private:
  AtomHandle first_,second_;
  Real len_, xi_, omega, xi_;
};


}} // ns

#endif
