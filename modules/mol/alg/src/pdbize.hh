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
#ifndef OST_MOL_ALG_PDBIZE_HH
#define OST_MOL_ALG_PDBIZE_HH

#include <ost/mol/entity_view.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/seq/sequence_list.hh>
#include "module_config.hh"

namespace ost { namespace mol { namespace alg {


extern const char* POLYPEPTIDE_CHAIN_NAMES;
extern const char* LIGAND_CHAIN_NAME;
extern const char* WATER_CHAIN_NAME;

class DLLEXPORT_OST_MOL_ALG PDBize {
public:
  explicit PDBize(int min_polymer_size=10):
    peptide_min_size_(min_polymer_size),
    nucleicacid_min_size_(min_polymer_size),
    saccharide_min_size_(min_polymer_size), ent_(mol::CreateEntity()),
    curr_chain_name_(POLYPEPTIDE_CHAIN_NAMES), needs_adjustment_(false),
    last_rnum_(0)
  {}

  explicit PDBize(int peptide_min_size,
                  int nucleicacid_min_size,
                  int saccharide_min_size):
    peptide_min_size_(peptide_min_size),
    nucleicacid_min_size_(nucleicacid_min_size),
    saccharide_min_size_(saccharide_min_size), ent_(mol::CreateEntity()),
    curr_chain_name_(POLYPEPTIDE_CHAIN_NAMES), needs_adjustment_(false),
    last_rnum_(0)
  {}

  void Add(mol::EntityView asu, const geom::Mat4List& transforms,
           seq::SequenceList seqres);

  EntityHandle Finish(bool shift_to_fit=true);
private:
  int                                   peptide_min_size_;
  int                                   nucleicacid_min_size_;
  int                                   saccharide_min_size_;
  EntityHandle                          ent_;
  ChainHandle                           ligand_chain_;
  ChainHandle                           water_chain_;
  const char*                           curr_chain_name_;
  bool                                  needs_adjustment_;
  int                                   last_rnum_;
  ResNum                                last_water_rnum_;
  std::map<ResidueHandle,ResidueHandle> dst_to_src_map_;
};

}}}
#endif

