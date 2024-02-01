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
#ifndef OST_IO_SAVE_ENTITY_HH
#define OST_IO_SAVE_ENTITY_HH

#include <ost/io/module_config.hh>

namespace ost { 
  
namespace mol {

class EntityHandle;
class EntityView;
}

namespace io {
  
/// \brief  Save entity handle
/// \param en is the entity to be saved
/// \param filename is the output filename
/// \param format
///         specifies what output format is used for export. If set to auto the
///         format is deduced from the file extension
void DLLEXPORT_OST_IO SaveEntity(const mol::EntityHandle& en, 
                                 const String& filename,
                                 const String& format="auto");
/// \brief  Save entity view
/// 
/// This method provides a convenient way to save only  a certain part of an 
/// entity, i.e. to filter out water and ligands or to create a backbone-only 
/// model.
/// 
/// \param en is the entity view to be saved.
/// \param filename is the output filename
/// \param format
///         specifies what output format is used for export. If set to auto the
///         format is deduced from the file extension
void DLLEXPORT_OST_IO SaveEntity(const mol::EntityView& en, 
                                const String& filename,
                                const String& format="auto");
}}

#endif
