//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2024 by the OpenStructure authors
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
#include <ost/conop/compound_lib.hh>

#include <ost/io/mol/mmcif_info.hh>
#include <ost/io/mol/io_profile.hh>
#include <ost/io/mol/star_writer.hh>

namespace ost { namespace io {


struct MMCifWriterEntity {

  MMCifWriterEntity() { }

  static MMCifWriterEntity FromPolymer(const String& entity_poly_type,
                                       const std::vector<String>& mon_ids,
                                       conop::CompoundLibPtr compound_lib);

  int GetAsymIdx(const String& asym_id) const;

  bool operator==(const MMCifWriterEntity& rhs) const {
       return (type == rhs.type)
           && (poly_type == rhs.poly_type)
           && (branch_type == rhs.branch_type)
           && (asym_ids == rhs.asym_ids)
           && (is_poly == rhs.is_poly)
           && (mon_ids == rhs.mon_ids)
           && (seq_olcs == rhs.seq_olcs)
           && (seq_can_olcs == rhs.seq_can_olcs)
           && (asym_alns == rhs.asym_alns);
  }

  bool operator!=(const MMCifWriterEntity& rhs) const {
    return !(*this == rhs);
  }

  // _entity.type
  String type;

  // _entity_poly.type
  String poly_type;  

  // __pdbx_entity_branch.type
  String branch_type;

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
  // irrelevant if is_poly is false. The assumption is that aligned residues
  // exactly match with the respective position in mon_ids.
  std::vector<std::vector<String> > asym_alns; 
};

class DLLEXPORT_OST_IO MMCifWriter : public StarWriter {
public:

  MMCifWriter(): structure_set_(false) { }

  virtual ~MMCifWriter() { }

  void SetStructure(const ost::mol::EntityHandle& ent, conop::CompoundLibPtr compound_lib,
                    bool mmcif_conform=true,
                    const std::vector<MMCifWriterEntity>& entity_info=std::vector<MMCifWriterEntity>());

  void SetStructure(const ost::mol::EntityView& ent, conop::CompoundLibPtr compound_lib,
                    bool mmcif_conform=true,
                    const std::vector<MMCifWriterEntity>& entity_info=std::vector<MMCifWriterEntity>());

  const std::vector<MMCifWriterEntity>& GetEntities() const { return entity_info_; }
  
private:

  void Setup();

  void Finalize(ost::conop::CompoundLibPtr compound_lib);

  std::vector<MMCifWriterEntity> entity_info_;
  StarWriterLoopPtr atom_type_;
  StarWriterLoopPtr atom_site_;
  StarWriterLoopPtr pdbx_poly_seq_scheme_;
  StarWriterLoopPtr entity_;
  StarWriterLoopPtr struct_asym_;
  StarWriterLoopPtr entity_poly_;
  StarWriterLoopPtr entity_poly_seq_;
  StarWriterLoopPtr chem_comp_;
  StarWriterLoopPtr pdbx_entity_branch_;
  bool structure_set_;
};

}} // ns

#endif
