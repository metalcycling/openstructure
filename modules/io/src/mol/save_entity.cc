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
#include <ost/log.hh>
#include <ost/io/io_manager.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/entity_view.hh>

#include "save_entity.hh"

namespace ost { namespace io {
  
void DLLEXPORT_OST_IO SaveEntity(const mol::EntityHandle& en, 
                                 const String& filename,
                                 const String& format) {
  SaveEntity(en.CreateFullView(), filename, format);
}

void DLLEXPORT_OST_IO SaveEntity(const mol::EntityView& en, 
                                const String& filename,
                                const String& format) {
  LOG_DEBUG("creating EntityIOHandle for " << filename);
  IOManager& manager=IOManager::Instance();
  EntityIOHandlerP ent_io=manager.FindEntityExportHandler(filename, format);

  if (!ent_io) {
    throw IOUnknownFormatException("No IO handler for file='"+filename+"', format='"+
                      format+"' found.");
  }
  LOG_DEBUG("calling export on entity io handle");

  ent_io->Export(en, filename);                                    
}

}}
