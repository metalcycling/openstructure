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
#include <sstream>
#include <ost/io/mol/pdb_str.hh>
#include <ost/io/mol/pdb_writer.hh>
#include <ost/io/mol/pdb_reader.hh>

namespace ost { namespace io {

String EntityToPDBString(const mol::EntityHandle& ent, const IOProfile& profile) {
  std::stringstream stream;
  PDBWriter writer(stream, profile);
  writer.Write(ent);
  return stream.str();
}

String EntityToPDBString(const mol::EntityView& ent, const IOProfile& profile) {
  std::stringstream stream;
  PDBWriter writer(stream, profile);
  writer.Write(ent);
  return stream.str();
}

mol::EntityHandle PDBStringToEntity(const String& pdb, const IOProfile& profile,
                                    bool process) {
  std::stringstream stream(pdb);
  PDBReader reader(stream, profile);
  mol::EntityHandle ent = mol::CreateEntity();
  reader.Import(ent);
  if(profile.processor && process) {
    profile.processor->Process(ent);
  }
  return ent;
}

}}
