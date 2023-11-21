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

#include <ost/mol/chem_class.hh>

#include <ost/io/mol/mmcif_writer.hh>



namespace {

  // generates as many chain names as you want (potentially multiple characters)
  struct ChainNameGenerator{
    ChainNameGenerator() { 
      chain_names = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz";
      n_chain_names = chain_names.size();
      indices.push_back(-1);
    }

    String Get() {
      int idx = indices.size() - 1;
      indices[idx] += 1;
      bool more_digits = false;
      while(idx >= 0) {
        if(indices[idx] >= n_chain_names) {
          indices[idx] = 0;
          if(idx>0) {
            indices[idx-1] += 1;
            --idx;
          } else {
            more_digits = true;
            break;
          }
        } else {
          break;
        }
      }
      if(more_digits) {
        indices.insert(indices.begin(), 0);
      }
      String ch_name(indices.size(), 'X');
      for(uint i = 0; i < indices.size(); ++i) {
        ch_name[i] = chain_names[indices[i]];
      }
      return ch_name;
    }

    void Reset() {
      indices.clear();
      indices.push_back(-1);
    }

    String chain_names;
    int n_chain_names;
    std::vector<int> indices;
  };
}

namespace ost { namespace io {

MMCifWriter::MMCifWriter(const String& filename, const IOProfile& profile):
  StarWriter(filename),
  profile_(profile),
  atom_site_(NULL) { }

MMCifWriter::~MMCifWriter() {
  if(atom_site_ != NULL) {
    delete atom_site_;
  }
}

void MMCifWriter::Process_atom_site(const ost::mol::EntityHandle& ent) {

  this->Setup_atom_site_();

  std::vector<std::vector<ost::mol::ResidueHandle> > L_chains; // L_PEPTIDE_LINKING
  std::vector<std::vector<ost::mol::ResidueHandle> > D_chains; // D_PEPTIDE_LINKING
  std::vector<std::vector<ost::mol::ResidueHandle> > P_chains; // PEPTIDE_LINKING  
  std::vector<std::vector<ost::mol::ResidueHandle> > R_chains; // RNA_LINKING
  std::vector<std::vector<ost::mol::ResidueHandle> > S_chains; // DNA_LINKING
  std::vector<std::vector<ost::mol::ResidueHandle> > X_chains; // L_SACCHARIDE
  std::vector<std::vector<ost::mol::ResidueHandle> > Y_chains; // D_SACCHARIDE
  std::vector<std::vector<ost::mol::ResidueHandle> > W_chains; // WATER
  std::vector<ost::mol::ResidueHandle> N_chains; // NON_POLYMER (1 res per chain)
  std::vector<ost::mol::ResidueHandle> U_chains; // UNKNOWN (1 res per chain)

  ost::mol::ChainHandleList chain_list = ent.GetChainList();
  for(auto ch: chain_list) {

    ost::mol::ResidueHandleList res_list = ch.GetResidueList();

    // we don't just go for chain type here...
    // just think of PDB entries that have a polypeptide, water and a ligand
    // in the same chain...

    bool has_l_peptide_linking = false;
    bool has_d_peptide_linking = false;
    bool has_peptide_linking = false;
    bool has_rna_linking = false;
    bool has_dna_linking = false;
    bool has_l_saccharide = false;
    bool has_d_saccharide = false;
    bool has_water = false;
    for(auto res: res_list) {

      // Peptide chains must not mix L_PEPTIDE_LINKING AND D_PEPTIDE_LINKING
      if(res.GetChemClass() == ost::mol::ChemClass::D_PEPTIDE_LINKING) {
        if(has_l_peptide_linking) {
          throw ost::io::IOException("Cannot write mmCIF when same chain "
                                     "contains D- and L-peptides");
        }
        has_d_peptide_linking = true;
      }
      if(res.GetChemClass() == ost::mol::ChemClass::L_PEPTIDE_LINKING) {
        if(has_d_peptide_linking) {
          throw ost::io::IOException("Cannot write mmCIF when same chain "
                                     "contains D- and L-peptides");
        }
        has_l_peptide_linking = true;
      }

      if(res.GetChemClass() == ost::mol::ChemClass::PEPTIDE_LINKING) {
        has_peptide_linking = true;
      }

      if(res.GetChemClass() == ost::mol::ChemClass::RNA_LINKING) {
        has_rna_linking = true;
      }

      if(res.GetChemClass() == ost::mol::ChemClass::DNA_LINKING) {
        has_dna_linking = true;
      }

      if(res.GetChemClass() == ost::mol::ChemClass::L_SACCHARIDE) {
        has_l_saccharide = true;
      }

      if(res.GetChemClass() == ost::mol::ChemClass::D_SACCHARIDE) {
        has_d_saccharide = true;
      }

      if(res.GetChemClass() == ost::mol::ChemClass::WATER) {
        has_water = true;
      }
    }

    // if there is any L-peptide or D-peptide, all peptides without
    // chiral center get assigned to it. No need for specific chain.
    if(has_l_peptide_linking || has_d_peptide_linking) {
      has_peptide_linking = false;
    }

    if(has_l_peptide_linking) {
      L_chains.push_back(ost::mol::ResidueHandleList());
    }

    if(has_d_peptide_linking) {
      D_chains.push_back(ost::mol::ResidueHandleList());
    }

    if(has_peptide_linking) {
      P_chains.push_back(ost::mol::ResidueHandleList());
    }

    if(has_rna_linking) {
      R_chains.push_back(ost::mol::ResidueHandleList());
    }

    if(has_dna_linking) {
      S_chains.push_back(ost::mol::ResidueHandleList());
    }

    if(has_l_saccharide) {
      X_chains.push_back(ost::mol::ResidueHandleList());
    }

    if(has_d_saccharide) {
      Y_chains.push_back(ost::mol::ResidueHandleList());
    }

    if(has_water) {
      W_chains.push_back(ost::mol::ResidueHandleList());
    }

    for(auto res: res_list) {
      if(res.GetChemClass().IsPeptideLinking()) {
        if(has_l_peptide_linking) {
          L_chains.back().push_back(res);
        } else if(has_d_peptide_linking) {
          D_chains.back().push_back(res);
        } else {
          P_chains.back().push_back(res);
        }
      } else if(res.GetChemClass() == ost::mol::ChemClass::RNA_LINKING) {
        R_chains.back().push_back(res);
      } else if(res.GetChemClass() == ost::mol::ChemClass::DNA_LINKING) {
        S_chains.back().push_back(res);
      } else if(res.GetChemClass() == ost::mol::ChemClass::L_SACCHARIDE) {
        X_chains.back().push_back(res);
      } else if(res.GetChemClass() == ost::mol::ChemClass::D_SACCHARIDE) {
        Y_chains.back().push_back(res);
      } else if(res.GetChemClass() == ost::mol::ChemClass::WATER) {
        W_chains.back().push_back(res);
      } else if(res.GetChemClass() == ost::mol::ChemClass::NON_POLYMER) {
        N_chains.push_back(res);
      } else if(res.GetChemClass() == ost::mol::ChemClass::UNKNOWN) {
        U_chains.push_back(res);
      } else {
        // TODO: make error message more insightful...
        throw ost::io::IOException("Unsupported chem class...");
      }
    }
  }

  ChainNameGenerator chain_name_gen;

  // process L_PEPTIDE_LINKING
  for(auto res_list: L_chains) {
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process D_PEPTIDE_LINKING
  for(auto res_list: D_chains) {
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process PEPTIDE_LINKING
  for(auto res_list: P_chains) {
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process RNA_LINKING
  for(auto res_list: R_chains) {
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process DNA_LINKING
  for(auto res_list: S_chains) {
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process L_SACHARIDE
  for(auto res_list: X_chains) {
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process D_SACHARIDE
  for(auto res_list: Y_chains) {
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process WATER
  for(auto res_list: W_chains) {
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process NON_POLYMER
  for(auto res: N_chains) {
    ost::mol::ResidueHandleList res_list;
    res_list.push_back(res);
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // process UNKNOWN
  for(auto res: N_chains) {
    ost::mol::ResidueHandleList res_list;
    res_list.push_back(res);
    String chain_name = chain_name_gen.Get();
    Feed_atom_site_(chain_name, 0, res_list);
  }

  // finalize
  this->Push(atom_site_);
}

void MMCifWriter::Setup_atom_site_() {
  StarLoopDesc desc;
  desc.SetCategory("_atom_site");
  desc.Add("group_PDB");
  desc.Add("type_symbol");
  desc.Add("label_atom_id");
  desc.Add("label_comp_id");
  desc.Add("label_asym_id");
  desc.Add("label_entity_id");
  desc.Add("label_seq_id");
  desc.Add("label_alt_id");
  desc.Add("Cartn_x");
  desc.Add("Cartn_y");
  desc.Add("Cartn_z");
  desc.Add("occupancy");
  desc.Add("B_iso_or_equiv");
  desc.Add("auth_seq_id");
  desc.Add("auth_asym_id");
  desc.Add("id");
  atom_site_ = new StarLoop(desc);
}

void MMCifWriter::Feed_atom_site_(const String& label_asym_id,
                                  int label_entity_id,
                                  const ost::mol::ResidueHandleList& res_list) {

  int label_seq_id = 1;
  for(auto res: res_list) {
    String comp_id = res.GetName();
    ost::mol::AtomHandleList at_list = res.GetAtomList();
    String auth_asym_id = res.GetChain().GetName();
    String auth_seq_id = res.GetNumber().AsString();
    for(auto at: at_list) {
      std::vector<StarLoopDataItemDO> at_data;
      // group_PDB
      if(at.IsHetAtom()) {
        at_data.push_back(StarLoopDataItemDO("HETATM"));
      } else {
        at_data.push_back(StarLoopDataItemDO("ATOM"));
      }
      // type_symbol
      at_data.push_back(StarLoopDataItemDO(at.GetElement()));
      // label_atom_id
      at_data.push_back(StarLoopDataItemDO(at.GetName()));
      // label_comp_id
      at_data.push_back(StarLoopDataItemDO(comp_id));
      // label_asym_id
      at_data.push_back(StarLoopDataItemDO(label_asym_id));
      // label_entity_id
      at_data.push_back(StarLoopDataItemDO(label_entity_id));
      // label_seq_id
      at_data.push_back(StarLoopDataItemDO(label_seq_id));
      // label_alt_id
      at_data.push_back(StarLoopDataItemDO("."));
      // Cartn_x
      at_data.push_back(StarLoopDataItemDO(at.GetPos().GetX(), 3));
      // Cartn_y
      at_data.push_back(StarLoopDataItemDO(at.GetPos().GetY(), 3));
      // Cartn_z
      at_data.push_back(StarLoopDataItemDO(at.GetPos().GetZ(), 3));
      // occupancy
      at_data.push_back(StarLoopDataItemDO(at.GetOccupancy(), 2));
      // B_iso_or_equiv
      at_data.push_back(StarLoopDataItemDO(at.GetBFactor(), 2));
      // auth_seq_id
      at_data.push_back(StarLoopDataItemDO(auth_seq_id));
      // auth_asym_id
      at_data.push_back(StarLoopDataItemDO(auth_asym_id));
      // id
      at_data.push_back(StarLoopDataItemDO(atom_site_->GetN()));
      atom_site_->AddData(at_data);
    }
    ++label_seq_id;
  }
}

}} // ns
