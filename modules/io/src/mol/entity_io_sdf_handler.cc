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
/*
  Author: Tobias Schmidt
 */

#include <ost/log.hh>
#include <ost/io/mol/sdf_writer.hh>
#include <ost/io/mol/sdf_reader.hh>

#include "entity_io_sdf_handler.hh"

namespace ost { namespace io {


bool EntityIOSDFHandler::RequiresProcessor() const
{
  return false;
}

void EntityIOSDFHandler::Import(mol::EntityHandle& ent,
                                std::istream& instream)
{
  SDFReader reader(instream);
  reader.Import(ent);
}

void EntityIOSDFHandler::Import(mol::EntityHandle& ent,
                                const boost::filesystem::path& loc)
{
  SDFReader reader(loc);
  reader.Import(ent);
}

void EntityIOSDFHandler::Export(const mol::EntityView& ent,
                                const boost::filesystem::path& loc) const {

  SDFWriter writer(loc);
  writer.Write(ent);
}

void EntityIOSDFHandler::Export(const mol::EntityView& ent, 
                    std::ostream& outstream) const
{
  SDFWriter writer(outstream);
  writer.Write(ent);
}

namespace {

bool sdf_handler_is_responsible_for(const boost::filesystem::path& loc,
                                    const String& type) {
  if(type=="auto") {
	String match_suf_string=loc.string();
    std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
    if(detail::FilenameEndsWith(match_suf_string,".sdf") || detail::FilenameEndsWith(match_suf_string,".sdf.gz")) {
      return true;
    }

  } else if(type=="sdf") {
    return true;
  }

  return false;
}

}

bool EntityIOSDFHandler::ProvidesImport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return sdf_handler_is_responsible_for(loc, type);
}

bool EntityIOSDFHandler::ProvidesExport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return sdf_handler_is_responsible_for(loc, type);
}

mol::EntityHandle LoadSDF(const String& file_name) {

    mol::EntityHandle ent_handle=mol::CreateEntity();
    EntityIOSDFHandler ent_io;
    ent_io.Import(ent_handle,file_name);

    return ent_handle;
}

void SaveSDF(const mol::EntityHandle& ent, const String& file_name) {
  SDFWriter writer(file_name);
  writer.Write(ent);
}

void SaveSDF(const mol::EntityView& ent, const String& file_name) {
  SDFWriter writer(file_name);
  writer.Write(ent);
}

}} // ns

