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
#ifndef OST_IO_PDB_STR_HH
#define OST_IO_PDB_STR_HH

#include <ost/io/module_config.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/io/mol/io_profile.hh>

namespace ost { namespace io {

// Saves an entity to a string in PDB format.
String DLLEXPORT_OST_IO
EntityToPDBString(const mol::EntityHandle& ent, const IOProfile& profile=IOProfile());

String DLLEXPORT_OST_IO
EntityToPDBString(const mol::EntityView& ent, const IOProfile& profile=IOProfile());

mol::EntityHandle DLLEXPORT_OST_IO
PDBStringToEntity(const String& pdb, const IOProfile& profile=IOProfile(),
                  bool process=false);

}}

#endif
