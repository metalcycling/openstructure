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
  Author: Marco Biasini, Ansgar Philippsen
 */
#include <iostream>
#include <sstream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <ost/log.hh>
#include <ost/io/io_exception.hh>
#include "entity_io_mmcif_handler.hh"
#include <ost/profile.hh>
#include <ost/io/mol/mmcif_reader.hh>
namespace ost { namespace io {

using boost::format;


bool EntityIOMMCIFHandler::RequiresProcessor() const
{
  return true;
}


void EntityIOMMCIFHandler::Export(const mol::EntityView& ent,
                                const boost::filesystem::path& loc) const 
{
}

void EntityIOMMCIFHandler::Import(mol::EntityHandle& ent, 
                                std::istream& stream)
{
  MMCifReader reader(stream,ent,  
                     IOProfileRegistry::Instance().GetDefault());
  reader.Parse();
  // This is a hack: the reader is not available in
  // file_loader.cc/ FileLoader::TryLoadEntity, but this is where the Conop
  // processor is called to establish covalent bonds. Since RequiresProcessor()
  // always returns true, the mechanism always tries to establish bonds, so we
  // do that here.
  MMCifInfo info = reader.GetInfo();
  info.ConnectBranchLinks();
}

void EntityIOMMCIFHandler::Export(const mol::EntityView& ent,
                                std::ostream& stream) const 
{
}

void EntityIOMMCIFHandler::Import(mol::EntityHandle& ent,
                                const boost::filesystem::path& loc)
{
  std::string filename=loc.string();
  MMCifReader reader(filename, ent,  
                     IOProfileRegistry::Instance().GetDefault());
  reader.Parse();
  // This is a hack: the reader is not available in
  // file_loader.cc/ FileLoader::TryLoadEntity, but this is where the Conop
  // processor is called to establish covalent bonds. Since RequiresProcessor()
  // always returns true, the mechanism always tries to establish bonds, so we
  // do that here.
  MMCifInfo info = reader.GetInfo();
  info.ConnectBranchLinks();
}

bool EntityIOMMCIFHandler::ProvidesImport(const boost::filesystem::path& loc,
                                        const String& type)
{
  if (type=="auto") {
    String match_suf_string=loc.string();
    std::transform(match_suf_string.begin(),match_suf_string.end(),
                   match_suf_string.begin(),tolower);
    if(detail::FilenameEndsWith(match_suf_string,".cif") ||
       detail::FilenameEndsWith(match_suf_string,".cif.gz")){
      return true;
    }
  }
  return type=="cif";
}

bool EntityIOMMCIFHandler::ProvidesExport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return false;
}


}} // ns


