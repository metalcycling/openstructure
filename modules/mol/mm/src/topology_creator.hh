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

#ifndef OST_MM_TOPOLOGY_CREATOR_HH
#define OST_MM_TOPOLOGY_CREATOR_HH

#include <vector>
#include <map>
#include <set>
#include <math.h>

#include <ost/message.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/mm/forcefield.hh>
#include <ost/mol/mm/buildingblock.hh>
#include <ost/mol/mm/block_modifiers.hh>
#include <ost/mol/mm/interaction.hh>
#include <ost/mol/mm/settings.hh>
#include <ost/mol/mm/modeller.hh>
#include <ost/mol/mm/index.hh>
#include <ost/mol/mm/topology.hh>


namespace ost { namespace mol{ namespace mm{

class TopologyCreator;
typedef boost::shared_ptr<TopologyCreator> TopologyCreatorPtr;

class TopologyCreator {
public:
  static TopologyPtr Create(ost::mol::EntityHandle& ent, 
                            const SettingsPtr settings);

};


}}}//ns

#endif
