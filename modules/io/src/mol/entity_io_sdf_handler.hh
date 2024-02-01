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
#ifndef OST_IO_ENTITY_IO_PLUGIN_SDF_H
#define OST_IO_ENTITY_IO_PLUGIN_SDF_H

#include <ost/io/mol/entity_io_handler.hh>

namespace ost { namespace io {

class DLLEXPORT_OST_IO EntityIOSDFHandler: public EntityIOHandler {
public:
  virtual void Import(mol::EntityHandle& ent, 
                      const boost::filesystem::path& loc);
  virtual void Export(const mol::EntityView& ent, 
                      const boost::filesystem::path& loc) const;
  // import data from provided stream
  virtual void Import(mol::EntityHandle& ent, 
                      std::istream& stream);
  // export data from entity view to provided stream
  virtual void Export(const mol::EntityView& ent, 
                      std::ostream& stream) const;
  static bool ProvidesImport(const boost::filesystem::path& loc, 
                             const String& format="auto");
  static bool ProvidesExport(const boost::filesystem::path& loc, 
                             const String& format="auto");
                             
  virtual bool RequiresProcessor() const;

  static String GetFormatName() { return String("Sdf"); }
  static String GetFormatDescription() { return String("Structure-data format from Molecular Design Limited"); }
};

typedef EntityIOHandlerFactory<EntityIOSDFHandler> EntityIOSDFHandlerFactory;

mol::EntityHandle DLLEXPORT_OST_IO LoadSDF(const String& file_name);

}} // ns

#endif
