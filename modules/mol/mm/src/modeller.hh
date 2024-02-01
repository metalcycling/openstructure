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

#ifndef OST_MM_MODELLER_HH
#define OST_MM_MODELLER_HH


#include <vector>
#include <set>

#include <ost/platform.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/geom/vec3.hh>
#include <ost/geom/composite.hh>
#include <ost/geom/transform.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/message.hh>
#include <ost/mol/bounding_box.hh>
#include <ost/mol/mm/buildingblock.hh>
#include <ost/mol/mm/forcefield.hh>
#include <ost/mol/mm/settings.hh>
#include <ost/mol/mm/index.hh>
#include <ost/mol/mm/topology.hh>
#include <ost/conop/heuristic.hh>
#include <ost/mol/spatial_organizer.hh>


namespace ost { namespace mol{ namespace mm{

class Modeller{
public:

  static void GenerateDisulfidBonds(ost::mol::EntityHandle& handle);

  //may sound pretty stupid...
  //But this is necessary if we want to compare our energies with
  //the energies calculated by gromacs. When ost writes and entity
  //down to disk, it uses 3 digits precision.
  //If we feed this into gromacs, it also produces a topology
  //with three digits precision. The problem is, that gromacs
  //uses nm, therefore the precision gets lowered tenfold.
  //To take care of that it is necessary to lower our precision...
  static void LowerPrecision(ost::mol::EntityHandle& handle);

  static void AssignPDBNaming(ost::mol::EntityHandle& handle);

  static void AssignGromacsNaming(ost::mol::EntityHandle& handle);
};


}}}

#endif
