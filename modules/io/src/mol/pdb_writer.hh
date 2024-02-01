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
#ifndef OST_IO_PDB_WRITER_HH
#define OST_IO_PDB_WRITER_HH


/*
  Author: Marco Biasini
 */
#include <fstream>

#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <ost/io/module_config.hh>
#include <ost/io/formatted_line.hh>
#include <ost/io/mol/io_profile.hh>

namespace ost { 
  
namespace mol {

class EntityView;
class EntityHandle;

}

namespace io {


class DLLEXPORT_OST_IO PDBWriter {
  typedef boost::iostreams::filtering_stream<boost::iostreams::output> OutStream;
public:
  PDBWriter(const String& filename,
            const IOProfile& profile);
  PDBWriter(const boost::filesystem::path& filename,
            const IOProfile& profile);
  PDBWriter(std::ostream& outstream, const IOProfile& profile);
  void SetWriteMultiModel(bool flag) { multi_model_=flag; }
  bool GetWriteMultiModel() const { return multi_model_; }
  void SetIsPQR(bool flag) { is_pqr_=flag; }
  bool IsPQR() const { return is_pqr_; }
  void Write(const mol::EntityView& ent);
  void Write(const mol::EntityHandle& ent);
  
  void Write(const mol::AtomHandleList& atoms);

  ~PDBWriter();
private:
  template <typename H>
  void WriteModel(H ent);
  
  void WriteModelLeader();
  void WriteModelTrailer();
  std::ofstream       outfile_;
  std::ostream&       outstream_;
  int                 mol_count_;
  std::map<long, int> atom_indices_;
  FormattedLine       line_;
  bool                multi_model_;
  bool                charmm_style_;
  bool                is_pqr_;
  IOProfile           profile_;
  String              filename_;
  OutStream           out_;
};
 
}}

#endif
