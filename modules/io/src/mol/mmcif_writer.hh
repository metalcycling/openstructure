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


struct MMCifWriterEntity {

  int GetAsymIdx(const String& asym_id) const;

  // _entity.type
  String type;

  // _entity_poly.type
  String poly_type;  

  // Names of chains in AU that are assigned to this entity
  std::vector<String> asym_ids;

  // in principle type == "polymer"
  bool is_poly;

  // SEQRES... kind of... internally we're not working on one letter codes
  // etc. but on full compound names. Only one element if is_poly is false.
  std::vector<String> mon_ids;

  // The respective strings for pdbx_seq_one_letter_code
  // Irrelevant if is_poly is false.
  std::vector<String> seq_olcs;

  // same for pdbx_seq_one_letter_code_can
  // Irrelevant if is_poly is false.
  std::vector<String> seq_can_olcs;

  // One alignment to mon_ids for each element in asym_ids, i.e. SEQRES-ATOMSEQ
  // alignment. Contains "-" for residues that are missing in ATOMSEQ.
  // irrelevant if is_poly is false.
  std::vector<std::vector<String> > asym_alns; 
};


class DLLEXPORT_OST_IO MMCifWriter : public StarWriter {
public:

  MMCifWriter(const String& filename, const IOProfile& profile);

  virtual ~MMCifWriter() { }

  void SetStructure(const ost::mol::EntityHandle& ent, bool mmcif_conform=true);

private:
  IOProfile profile_;
  std::vector<MMCifWriterEntity> entity_info_;
  StarWriterLoopPtr atom_type_;
  StarWriterLoopPtr atom_site_;
  StarWriterLoopPtr pdbx_poly_seq_scheme_;
  StarWriterLoopPtr entity_;
  StarWriterLoopPtr struct_asym_;
  StarWriterLoopPtr entity_poly_;
  StarWriterLoopPtr entity_poly_seq_;
  StarWriterLoopPtr chem_comp_;
};

}} // ns

#endif
