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

  // internal object with all info to fill chem_comp_ category
  struct CompInfo {
    String type;
  };

  inline String chem_class_to_chem_comp_type(char chem_class) {
    String type = "";
    switch(chem_class) {
      case 'P': {
        type = "peptide linking";
        break;
      }
      case 'D': {
        type = "D-peptide linking";
        break;
      }
      case 'L': {
        type = "L-peptide linking";
        break;
      }
      case 'R': {
        type = "RNA linking";
        break;
      }
      case 'S': {
        type = "DNA linking";
        break;
      }
      case 'N': {
        type = "non-polymer";
        break;
      }
      case 'X': {
        type = "L-saccharide";
        break;
      }
      case 'Y': {
        type = "D-saccharide";
        break;
      }
      case 'Z': {
        type = "saccharide";
        break;
      }
      case 'W': {
        type = "non-polymer";
        break;
      }
      case 'U': {
        type = "non-polymer";
        break;
      }
      default: {
        std::stringstream err;
        err << "Invalid chem class: "<<chem_class;
        throw ost::io::IOException(err.str());

      }
    }
    return type;
  }

  inline String chem_class_to_entity_poly_type(char chem_class) {
    String type = "";
    switch(chem_class) {
      case 'P': {
        type = "polypeptide(L)";
        break;
      }
      case 'D': {
        type = "polypeptide(D)";
        break;
      }
      case 'L': {
        type = "polypeptide(L)";
        break;
      }
      case 'R': {
        type = "polyribonucleotide";
        break;
      }
      case 'S': {
        type = "polydeoxyribonucleotide";
        break;
      }
      case 'X': {
        type = "polysaccharide(L)";
        break;
      }
      case 'Y': {
        type = "polysaccharide(D)";
        break;
      }
      default: {
        type = "other";
      }
    }
    return type;
  }

  void Setup_chem_comp_(const ost::mol::ResidueHandleList& res_list,
                        std::map<String, CompInfo>& comp_infos) {
    for(auto res: res_list) {
      String res_name = res.GetName();
      String type = chem_class_to_chem_comp_type(res.GetChemClass());
      auto it = comp_infos.find(res_name);
      if(it != comp_infos.end()) {
        // check whether type is consistent
        if(it->second.type != type) {
          throw ost::io::IOException("There can be only one");
        }
      } else {
        CompInfo info;
        info.type = type;
        comp_infos[res_name] = info;
      }
    }
  }

  // internal object with all info to fill entity_, struct_asym_,
  // entity_poly_seq_ categories
  struct EntityInfo {
    char chem_class; // all residues of this entity have this ChemClass
    String poly_type; // relevant for _entity_poly
    std::vector<String> asym_ids; // relevant for _struct_asym.id
    std::vector<String> mon_ids; // relevant for _entity_poly_seq.mon_id
    bool is_poly; // in principle mon_ids.size() > 1
  };

  int Setup_entity_(const String& asym_chain_name,
                    char chem_class,
                    const ost::mol::ResidueHandleList& res_list,
                    std::vector<EntityInfo>& entity_infos) {


    // deal with water
    if(chem_class == ost::mol::ChemClass::WATER) {
      for(size_t i = 0; i < entity_infos.size(); ++i) {
        if(entity_infos[i].chem_class == ost::mol::ChemClass::WATER) {
          entity_infos[i].asym_ids.push_back(asym_chain_name);
          return i;
        }
      }
      int entity_idx = entity_infos.size();
      entity_infos.push_back(EntityInfo());
      entity_infos.back().chem_class = ost::mol::ChemClass::WATER;
      entity_infos.back().asym_ids.push_back(asym_chain_name);
      entity_infos.back().mon_ids.push_back("HOH");
      entity_infos.back().poly_type = "";
      entity_infos.back().is_poly = false;
      return entity_idx; 
    }

    // deal with NON_POLYMER
    if(chem_class == ost::mol::ChemClass::NON_POLYMER) {
      for(size_t i = 0; i < entity_infos.size(); ++i) {
        if(entity_infos[i].chem_class == ost::mol::ChemClass::NON_POLYMER &&
           res_list[0].GetName() == entity_infos[i].mon_ids[0]) {
          entity_infos[i].asym_ids.push_back(asym_chain_name);
          return i;
        }
      }
      int entity_idx = entity_infos.size();
      entity_infos.push_back(EntityInfo());
      entity_infos.back().chem_class = ost::mol::ChemClass::NON_POLYMER;
      entity_infos.back().asym_ids.push_back(asym_chain_name);
      entity_infos.back().mon_ids.push_back(res_list[0].GetName());
      entity_infos.back().poly_type = "";
      entity_infos.back().is_poly = false;
      return entity_idx;
    }

    // deal with UNKNOWN
    if(chem_class == ost::mol::ChemClass::UNKNOWN) {
      for(size_t i = 0; i < entity_infos.size(); ++i) {
        if(entity_infos[i].chem_class == ost::mol::ChemClass::UNKNOWN &&
           res_list[0].GetName() == entity_infos[i].mon_ids[0]) {
          entity_infos[i].asym_ids.push_back(asym_chain_name);
          return i;
        }
      }
      int entity_idx = entity_infos.size();
      entity_infos.push_back(EntityInfo());
      entity_infos.back().chem_class = ost::mol::ChemClass::UNKNOWN;
      entity_infos.back().asym_ids.push_back(asym_chain_name);
      entity_infos.back().mon_ids.push_back(res_list[0].GetName());
      entity_infos.back().poly_type = "";
      entity_infos.back().is_poly = false;
      return entity_idx;
    }

    // with the current code, the following chem classes are considered
    // polymers: PEPTIDE_LINKING, D_PEPTIDE_LINKING, L_PEPTIDE_LINKING
    // RNA_LINKING, DNA_LINKING, L_SACCHARIDE, D_SACCHARIDE, SACCHARIDE
    // They're also considered polymers even if only one residue is there
    // Needs checking...

    std::vector<String> mon_ids;
    for(auto res : res_list) {
      mon_ids.push_back(res.GetName());
    }

    // check whether we already have that entity
    // right now we're just looking for exact matches in chem_class and
    // mon_ids (i.e. sequence). 
    int entity_idx = -1;
    for(size_t i = 0; i < entity_infos.size(); ++i) {
      if(entity_infos[i].chem_class == chem_class &&
         entity_infos[i].mon_ids == mon_ids) {
        entity_idx = i;
        break;
      }
    }

    if(entity_idx != -1) {
      entity_infos[entity_idx].asym_ids.push_back(asym_chain_name);
    } else {
      entity_idx = entity_infos.size();
      entity_infos.push_back(EntityInfo());
      entity_infos.back().chem_class = chem_class;
      entity_infos.back().asym_ids.push_back(asym_chain_name);
      entity_infos.back().mon_ids = mon_ids;
      entity_infos.back().poly_type = chem_class_to_entity_poly_type(chem_class);
      entity_infos.back().is_poly = entity_infos.back().mon_ids.size() > 1;
    }
    return entity_idx;
  }

  ost::io::StarLoop* Setup_atom_type_ptr() {
    ost::io::StarLoopDesc desc;
    desc.SetCategory("_atom_type");
    desc.Add("symbol");
    ost::io::StarLoop* sl = new ost::io::StarLoop(desc);
    return sl;
  }

  ost::io::StarLoop* Setup_atom_site_ptr() {
    ost::io::StarLoopDesc desc;
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
    desc.Add("pdbx_PDB_ins_code");
    ost::io::StarLoop* sl = new ost::io::StarLoop(desc);
    return sl;
  }

  ost::io::StarLoop* Setup_pdbx_poly_seq_scheme_ptr() {
    ost::io::StarLoopDesc desc;
    desc.SetCategory("_pdbx_poly_seq_scheme");
    desc.Add("asym_id");
    desc.Add("entity_id");
    desc.Add("mon_id");
    desc.Add("seq_id");
    desc.Add("pdb_strand_id");
    desc.Add("pdb_seq_num");
    desc.Add("pdb_ins_code");
    ost::io::StarLoop* sl = new ost::io::StarLoop(desc);
    return sl;
  }

  ost::io::StarLoop* Setup_entity_ptr() {
    ost::io::StarLoopDesc desc;
    desc.SetCategory("_entity");
    desc.Add("id");
    desc.Add("type");
    ost::io::StarLoop* sl = new ost::io::StarLoop(desc);
    return sl;
  }

  ost::io::StarLoop* Setup_struct_asym_ptr() {
    ost::io::StarLoopDesc desc;
    desc.SetCategory("_struct_asym");
    desc.Add("id");
    desc.Add("entity_id");
    ost::io::StarLoop* sl = new ost::io::StarLoop(desc);
    return sl;    
  }

  ost::io::StarLoop* Setup_entity_poly_ptr() {
    ost::io::StarLoopDesc desc;
    desc.SetCategory("_entity_poly");
    desc.Add("entity_id");
    desc.Add("type");
    ost::io::StarLoop* sl = new ost::io::StarLoop(desc);
    return sl;    
  }

  ost::io::StarLoop* Setup_entity_poly_seq_ptr() {
    ost::io::StarLoopDesc desc;
    desc.SetCategory("_entity_poly_seq");
    desc.Add("entity_id");
    desc.Add("mon_id");
    desc.Add("num");
    ost::io::StarLoop* sl = new ost::io::StarLoop(desc);
    return sl;    
  }

  ost::io::StarLoop* Setup_chem_comp_ptr() {
    ost::io::StarLoopDesc desc;
    desc.SetCategory("_chem_comp");
    desc.Add("id");
    desc.Add("type");
    ost::io::StarLoop* sl = new ost::io::StarLoop(desc);
    return sl;    
  }

  void Feed_atom_type_(ost::io::StarLoop* atom_type_ptr,
                       ost::io::StarLoop* atom_site_ptr) {
    // we're just extracting every type_symbol that we observed
    // in atom_site (this is a bit of circular stupidity...)
    std::set<String> symbols;
    int desc_size = atom_site_ptr->GetDesc().GetSize();
    int type_symbol_idx = atom_site_ptr->GetDesc().GetIndex("type_symbol");
    int N = atom_site_ptr->GetN();
    const std::vector<ost::io::StarLoopDataItemDO>& data = atom_site_ptr->GetData();
    for(int i = 0; i < N; ++i) {
      symbols.insert(data[i*desc_size + type_symbol_idx].GetValue());
    }
    std::vector<ost::io::StarLoopDataItemDO> atom_type_data;
    atom_type_data.push_back(ost::io::StarLoopDataItemDO(""));
    for(auto symbol: symbols) {
      atom_type_data[0] = ost::io::StarLoopDataItemDO(symbol);
      atom_type_ptr->AddData(atom_type_data);
    }
  }

  void Feed_pdbx_poly_seq_scheme(ost::io::StarLoop* pdbx_poly_seq_scheme_ptr,
                                 const String& label_asym_id,
                                 int label_entity_id,
                                 const ost::mol::ResidueHandleList& res_list) {

    std::vector<ost::io::StarLoopDataItemDO> data;
    data.push_back(ost::io::StarLoopDataItemDO(label_asym_id));
    data.push_back(ost::io::StarLoopDataItemDO(label_entity_id));
    data.push_back(ost::io::StarLoopDataItemDO(""));
    data.push_back(ost::io::StarLoopDataItemDO(0));
    data.push_back(ost::io::StarLoopDataItemDO(""));
    data.push_back(ost::io::StarLoopDataItemDO(0));
    data.push_back(ost::io::StarLoopDataItemDO(""));
    int label_seq_id = 1;
    for(auto res: res_list) {
      data[2] = ost::io::StarLoopDataItemDO(res.GetName());
      data[3] = ost::io::StarLoopDataItemDO(label_seq_id);
      data[4] = ost::io::StarLoopDataItemDO(res.GetChain().GetName());
      data[5] = ost::io::StarLoopDataItemDO(res.GetNumber().GetNum());
      char ins_code = res.GetNumber().GetInsCode();      
      if(ins_code == '\0') {
        data[6] = ost::io::StarLoopDataItemDO("");
      } else {
        String tmp = " ";
        tmp[0] = ins_code;
        data[6] = ost::io::StarLoopDataItemDO(tmp);
      }      
      pdbx_poly_seq_scheme_ptr->AddData(data);
      label_seq_id += 1;
    }
  }

  void Feed_atom_site_(ost::io::StarLoop* atom_site_ptr,
                       const String& label_asym_id,
                       int label_entity_id,
                       const ost::mol::ResidueHandleList& res_list,
                       bool is_poly) {
    int label_seq_id = 1;
    for(auto res: res_list) {
      String comp_id = res.GetName();
      ost::mol::AtomHandleList at_list = res.GetAtomList();
      String auth_asym_id = res.GetChain().GetName();
      String auth_seq_id = res.GetNumber().AsString();
      for(auto at: at_list) {
        std::vector<ost::io::StarLoopDataItemDO> at_data;
        // group_PDB
        if(at.IsHetAtom()) {
          at_data.push_back(ost::io::StarLoopDataItemDO("HETATM"));
        } else {
          at_data.push_back(ost::io::StarLoopDataItemDO("ATOM"));
        }
        // type_symbol
        at_data.push_back(ost::io::StarLoopDataItemDO(at.GetElement()));
        // label_atom_id
        at_data.push_back(ost::io::StarLoopDataItemDO(at.GetName()));
        // label_comp_id
        at_data.push_back(ost::io::StarLoopDataItemDO(comp_id));
        // label_asym_id
        at_data.push_back(ost::io::StarLoopDataItemDO(label_asym_id));
        // label_entity_id
        at_data.push_back(ost::io::StarLoopDataItemDO(label_entity_id));
        // label_seq_id
        if(is_poly) {
          at_data.push_back(ost::io::StarLoopDataItemDO(label_seq_id));
        } else {
          at_data.push_back(ost::io::StarLoopDataItemDO("."));
        }
        // label_alt_id
        at_data.push_back(ost::io::StarLoopDataItemDO("."));
        // Cartn_x
        at_data.push_back(ost::io::StarLoopDataItemDO(at.GetPos().GetX(), 3));
        // Cartn_y
        at_data.push_back(ost::io::StarLoopDataItemDO(at.GetPos().GetY(), 3));
        // Cartn_z
        at_data.push_back(ost::io::StarLoopDataItemDO(at.GetPos().GetZ(), 3));
        // occupancy
        at_data.push_back(ost::io::StarLoopDataItemDO(at.GetOccupancy(), 2));
        // B_iso_or_equiv
        at_data.push_back(ost::io::StarLoopDataItemDO(at.GetBFactor(), 2));
        // auth_seq_id
        at_data.push_back(ost::io::StarLoopDataItemDO(auth_seq_id));
        // auth_asym_id
        at_data.push_back(ost::io::StarLoopDataItemDO(auth_asym_id));
        // id
        at_data.push_back(ost::io::StarLoopDataItemDO(atom_site_ptr->GetN()));
        // pdbx_PDB_ins_code
        at_data.push_back(ost::io::StarLoopDataItemDO(""));
        atom_site_ptr->AddData(at_data);
      }
      ++label_seq_id;
    }
  }

  void Feed_entity_(ost::io::StarLoop* entity_ptr,
                    const std::vector<EntityInfo>& entity_info) {
    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      std::vector<ost::io::StarLoopDataItemDO> ent_data;
      // id
      ent_data.push_back(ost::io::StarLoopDataItemDO(entity_idx));
      // type
      ost::mol::ChemClass chem_class(entity_info[entity_idx].chem_class);
      if(chem_class.IsPeptideLinking() || chem_class.IsNucleotideLinking()) {
        ent_data.push_back(ost::io::StarLoopDataItemDO("polymer"));
      } else if(chem_class.IsWater()) {
        ent_data.push_back(ost::io::StarLoopDataItemDO("water"));
      } else if(chem_class == ost::mol::ChemClass::NON_POLYMER) {
        ent_data.push_back(ost::io::StarLoopDataItemDO("non-polymer"));        
      } else if(chem_class.IsSaccharide()) {
        // NOT SURE WHETHER THIS MAKES ANY SENSE!
        ent_data.push_back(ost::io::StarLoopDataItemDO("branched"));
      } else if(chem_class == ost::mol::ChemClass::UNKNOWN) {
        // NOT SURE WHETHER THIS MAKES ANY SENSE!
        ent_data.push_back(ost::io::StarLoopDataItemDO("non-polymer"));
      } else {
        throw ost::io::IOException("Entity type issue");
      }
      entity_ptr->AddData(ent_data);
    }
  }

  void Feed_struct_asym_(ost::io::StarLoop* struct_asym_ptr,
                         const std::vector<EntityInfo>& entity_info) {
    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      for(auto asym_id : entity_info[entity_idx].asym_ids) {
        std::vector<ost::io::StarLoopDataItemDO> asym_data;
        asym_data.push_back(ost::io::StarLoopDataItemDO(asym_id));
        asym_data.push_back(ost::io::StarLoopDataItemDO(entity_idx));
        struct_asym_ptr->AddData(asym_data);
      }
    }
  }

  void Feed_entity_poly_seq_(ost::io::StarLoop* entity_poly_seq_ptr,
                             const std::vector<EntityInfo>& entity_info) {
    // reuse data vector for efficiency
    std::vector<ost::io::StarLoopDataItemDO> entity_poly_seq_data;
    entity_poly_seq_data.push_back(ost::io::StarLoopDataItemDO(0));
    entity_poly_seq_data.push_back(ost::io::StarLoopDataItemDO("ALA"));
    entity_poly_seq_data.push_back(ost::io::StarLoopDataItemDO(1));

    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      if(entity_info[entity_idx].is_poly) {
        const std::vector<String>& mon_ids = entity_info[entity_idx].mon_ids;
        for(size_t mon_idx = 0; mon_idx < mon_ids.size(); ++mon_idx) {
          entity_poly_seq_data[0] = ost::io::StarLoopDataItemDO(entity_idx);
          entity_poly_seq_data[1] = ost::io::StarLoopDataItemDO(mon_ids[mon_idx]);
          entity_poly_seq_data[2] = ost::io::StarLoopDataItemDO(mon_idx+1);
          entity_poly_seq_ptr->AddData(entity_poly_seq_data);
        }
      }
    }
  }

  void Feed_entity_poly_(ost::io::StarLoop* entity_poly_ptr,
                         const std::vector<EntityInfo>& entity_info) {
    // reuse data vector for efficiency
    std::vector<ost::io::StarLoopDataItemDO> entity_poly_data;
    entity_poly_data.push_back(ost::io::StarLoopDataItemDO(0));
    entity_poly_data.push_back(ost::io::StarLoopDataItemDO("other"));
    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      if(entity_info[entity_idx].is_poly) {
        entity_poly_data[0] = ost::io::StarLoopDataItemDO(entity_idx);
        entity_poly_data[1] = ost::io::StarLoopDataItemDO(entity_info[entity_idx].poly_type);
        entity_poly_ptr->AddData(entity_poly_data);
      }
    }
  }

  void Feed_chem_comp_(ost::io::StarLoop* chem_comp_ptr,
                       const std::map<String, CompInfo>& comp_infos) {
    std::vector<ost::io::StarLoopDataItemDO> comp_data;
    comp_data.push_back(ost::io::StarLoopDataItemDO("ALA"));
    comp_data.push_back(ost::io::StarLoopDataItemDO("L-peptide linking"));
    for(auto it = comp_infos.begin(); it != comp_infos.end(); ++it) {
      comp_data[0] = it->first;
      comp_data[1] = it->second.type;
      chem_comp_ptr->AddData(comp_data);
    }
  }
}

namespace ost { namespace io {

MMCifWriter::MMCifWriter(const String& filename, const IOProfile& profile):
  StarWriter(filename),
  profile_(profile),
  atom_type_(NULL),
  atom_site_(NULL),
  pdbx_poly_seq_scheme_(NULL),
  entity_(NULL),
  struct_asym_(NULL),
  entity_poly_(NULL),
  entity_poly_seq_(NULL),
  chem_comp_(NULL) { }

MMCifWriter::~MMCifWriter() {
  if(atom_type_ != NULL) {
    delete atom_type_;
  }
  if(atom_site_ != NULL) {
    delete atom_site_;
  }
  if(pdbx_poly_seq_scheme_ != NULL) {
    delete pdbx_poly_seq_scheme_;
  }
  if(entity_ != NULL) {
    delete entity_;
  }
  if(struct_asym_ != NULL) {
    delete struct_asym_;
  }
  if(entity_poly_ != NULL) {
    delete entity_poly_;
  }
  if(entity_poly_seq_ != NULL) {
    delete entity_poly_seq_;
  }
  if(chem_comp_ != NULL) {
    delete chem_comp_;
  }
}

void MMCifWriter::SetEntity(const ost::mol::EntityHandle& ent) {

  // tabula rasa
  if(atom_type_ != NULL) {
    delete atom_type_;
  }
  if(atom_site_ != NULL) {
    delete atom_site_;
  }
  if(pdbx_poly_seq_scheme_ != NULL) {
    delete pdbx_poly_seq_scheme_;
  }
  if(entity_ != NULL) {
    delete entity_;
  }
  if(struct_asym_ != NULL) {
    delete struct_asym_;
  }
  if(entity_poly_ != NULL) {
    delete entity_poly_;
  }
  if(entity_poly_seq_ != NULL) {
    delete entity_poly_seq_;
  }
  if(chem_comp_ != NULL) {
    delete chem_comp_;
  }

  atom_type_ = Setup_atom_type_ptr();
  atom_site_ = Setup_atom_site_ptr();
  pdbx_poly_seq_scheme_ = Setup_pdbx_poly_seq_scheme_ptr();
  entity_ = Setup_entity_ptr();
  struct_asym_ = Setup_struct_asym_ptr();
  entity_poly_ = Setup_entity_poly_ptr();
  entity_poly_seq_ = Setup_entity_poly_seq_ptr();
  chem_comp_ = Setup_chem_comp_ptr();

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
  std::map<String, CompInfo> comp_infos;

  ost::mol::ChainHandleList chain_list = ent.GetChainList();
  for(auto ch: chain_list) {

    ost::mol::ResidueHandleList res_list = ch.GetResidueList();

    Setup_chem_comp_(res_list, comp_infos);

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
  std::vector<EntityInfo> entity_info;

  // process L_PEPTIDE_LINKING
  for(auto res_list: L_chains) {
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::L_PEPTIDE_LINKING,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process D_PEPTIDE_LINKING
  for(auto res_list: D_chains) {
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::D_PEPTIDE_LINKING,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process PEPTIDE_LINKING
  for(auto res_list: P_chains) {
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::PEPTIDE_LINKING,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process RNA_LINKING
  for(auto res_list: R_chains) {
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::RNA_LINKING,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process DNA_LINKING
  for(auto res_list: S_chains) {
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::DNA_LINKING,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process L_SACHARIDE
  for(auto res_list: X_chains) {
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::L_SACCHARIDE,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process D_SACHARIDE
  for(auto res_list: Y_chains) {
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::D_SACCHARIDE,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process WATER
  for(auto res_list: W_chains) {
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::WATER,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process NON_POLYMER
  for(auto res: N_chains) {
    ost::mol::ResidueHandleList res_list;
    res_list.push_back(res);
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::NON_POLYMER,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  // process UNKNOWN
  for(auto res: U_chains) {
    ost::mol::ResidueHandleList res_list;
    res_list.push_back(res);
    String chain_name = chain_name_gen.Get();
    int entity_id = Setup_entity_(chain_name,
                                  ost::mol::ChemClass::UNKNOWN,
                                  res_list,
                                  entity_info);
    Feed_atom_site_(atom_site_, chain_name, entity_id, res_list,
                    entity_info[entity_id].is_poly);
    if(entity_info[entity_id].is_poly) {
      Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme_, chain_name,
                                entity_id, res_list);
    }
  }

  Feed_entity_(entity_, entity_info);
  Feed_struct_asym_(struct_asym_, entity_info);
  Feed_entity_poly_(entity_poly_, entity_info);
  Feed_entity_poly_seq_(entity_poly_seq_, entity_info);
  Feed_chem_comp_(chem_comp_, comp_infos);
  Feed_atom_type_(atom_type_, atom_site_); 

  // finalize
  this->Push(chem_comp_);
  this->Push(entity_);
  this->Push(struct_asym_);
  this->Push(entity_poly_);
  this->Push(entity_poly_seq_);
  this->Push(pdbx_poly_seq_scheme_);
  this->Push(atom_type_);
  this->Push(atom_site_);
}

}} // ns
