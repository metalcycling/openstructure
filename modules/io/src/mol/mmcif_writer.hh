//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2023 by the OpenStructure authors
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
#ifndef OST_IO_MMCIF_WRITER_HH
#define OST_IO_MMCIF_WRITER_HH

#include <fstream>

#include <ost/mol/entity_handle.hh>

#include <ost/io/mol/mmcif_info.hh>
#include <ost/io/mol/io_profile.hh>
#include <ost/io/mol/star_writer.hh>

namespace ost { namespace io {

class DLLEXPORT_OST_IO MMCifWriter : public StarWriter {
public:

  MMCifWriter(const String& filename, const IOProfile& profile);

  virtual ~MMCifWriter();

  // initializes
  // atom_site_
  // entity_
  // struct_asym_
  // entity_poly_seq_
  void Process_atom_site(const ost::mol::EntityHandle& ent);

private:

  IOProfile profile_;
  StarLoop* atom_site_;
  StarLoop* entity_;
  StarLoop* struct_asym_;
  StarLoop* entity_poly_seq_;
};

 
}} // ns

#endif