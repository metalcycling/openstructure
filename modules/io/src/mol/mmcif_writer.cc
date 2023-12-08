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

  String GuessEntityPolyType(const ost::mol::ResidueHandleList& res_list) {

    // guesses _entity_poly.type based on residue chem classes

    // allowed values according to mmcif_pdbx_v50.dic:
    // - cyclic-pseudo-peptide 	
    // - other 	
    // - peptide nucleic acid 	
    // - polydeoxyribonucleotide 	
    // - polydeoxyribonucleotide/polyribonucleotide hybrid 	
    // - polypeptide(D) 	
    // - polypeptide(L) 	
    // - polyribonucleotide

    // this function won't identify cyclic-pseudo-peptides

    std::set<char> chem_classes;
    for(auto res: res_list) {
      chem_classes.insert(res.GetChemClass());
    }

    // check for polypeptide(L)
    if(chem_classes.size() == 1 &&
       chem_classes.find(ost::mol::ChemClass::L_PEPTIDE_LINKING) != chem_classes.end()) {
        return "polypeptide(L)";
    }

    if(chem_classes.size() == 2 &&
       chem_classes.find(ost::mol::ChemClass::L_PEPTIDE_LINKING) != chem_classes.end() &&
       chem_classes.find(ost::mol::ChemClass::PEPTIDE_LINKING) != chem_classes.end()) {
        return "polypeptide(L)";
    }

    // check for polypeptide(D)
    if(chem_classes.size() == 1 &&
       chem_classes.find(ost::mol::ChemClass::D_PEPTIDE_LINKING) != chem_classes.end()) {
        return "polypeptide(D)";
    }

    if(chem_classes.size() == 2 &&
       chem_classes.find(ost::mol::ChemClass::D_PEPTIDE_LINKING) != chem_classes.end() &&
       chem_classes.find(ost::mol::ChemClass::PEPTIDE_LINKING) != chem_classes.end()) {
        return "polypeptide(D)";
    }

    // check for polydeoxyribonucleotide
    if(chem_classes.size() == 1 &&
       chem_classes.find(ost::mol::ChemClass::DNA_LINKING) != chem_classes.end()) {
      return "polydeoxyribonucleotide";
    }

    // check for polyribonucleotide
    if(chem_classes.size() == 1 &&
       chem_classes.find(ost::mol::ChemClass::RNA_LINKING) != chem_classes.end()) {
      return "polyribonucleotide";
    }

    // check for polydeoxyribonucleotide/polyribonucleotide hybrid
    if(chem_classes.size() == 2 &&
       chem_classes.find(ost::mol::ChemClass::RNA_LINKING) != chem_classes.end() &&
       chem_classes.find(ost::mol::ChemClass::DNA_LINKING) != chem_classes.end()) {
      return "polydeoxyribonucleotide/polyribonucleotide hybrid";
    }

    // check for peptide nucleic acid
    bool peptide_linking = chem_classes.find(ost::mol::ChemClass::L_PEPTIDE_LINKING) != chem_classes.end() ||
                           chem_classes.find(ost::mol::ChemClass::D_PEPTIDE_LINKING) != chem_classes.end() ||
                           chem_classes.find(ost::mol::ChemClass::PEPTIDE_LINKING) != chem_classes.end();
    bool nucleotide_linking = chem_classes.find(ost::mol::ChemClass::DNA_LINKING) != chem_classes.end() ||
                              chem_classes.find(ost::mol::ChemClass::RNA_LINKING) != chem_classes.end();
    std::set<char> pepnuc_set;
    pepnuc_set.insert(ost::mol::ChemClass::L_PEPTIDE_LINKING);
    pepnuc_set.insert(ost::mol::ChemClass::D_PEPTIDE_LINKING);
    pepnuc_set.insert(ost::mol::ChemClass::PEPTIDE_LINKING);
    pepnuc_set.insert(ost::mol::ChemClass::DNA_LINKING);
    pepnuc_set.insert(ost::mol::ChemClass::RNA_LINKING);
    pepnuc_set.insert(chem_classes.begin(), chem_classes.end());
    if(peptide_linking && nucleotide_linking && pepnuc_set.size() == 5) {
      return "peptide nucleic acid";
    }

    return "other";
  }

  String GuessEntityType(const ost::mol::ResidueHandleList& res_list) {

    // guesses _entity.type based on residue chem classes

    // allowed values according to mmcif_pdbx_v50.dic:
    // - branched
    // - macrolide
    // - non-polymer
    // - polymer
    // - water

    // this function won't identify macrolid

    std::set<char> chem_classes;
    for(auto res: res_list) {
      chem_classes.insert(res.GetChemClass());
    }

    // check for water
    if(chem_classes.size() == 1 &&
       chem_classes.find(ost::mol::ChemClass::WATER) != chem_classes.end()) {
      return "water";
    }

    // check for non-polymer
    if(res_list.size() == 1) {
      return "non-polymer";
    }

    // check for branched
    std::set<char> sweet_set;
    sweet_set.insert(ost::mol::ChemClass::L_SACCHARIDE);
    sweet_set.insert(ost::mol::ChemClass::D_SACCHARIDE);
    sweet_set.insert(ost::mol::ChemClass::SACCHARIDE);
    // if the union of chem_classes and sweet_set has 3 elements, chem_classes
    // only has sugars.
    sweet_set.insert(chem_classes.begin(), chem_classes.end());
    if(sweet_set.size() == 3) {
      return "branched";
    }

    // DISCUSS THIS OVER A BEER...
    // when arriving here, we excluded the possibility of branched and single
    // residue chains.
    // BUT: entities must have at least 3 residues to be considered polymers
    // for now, we just set entities with 2 residues as non-polymer
    if(res_list.size() == 2) {
      return "non-polymer";
    }

    // If res_list represents a valid mmcif chain, chem_classes should only
    // contain peptide- and nucleotide linking items 
    for(auto it: chem_classes) {
      ost::mol::ChemClass chem_class(it);
      if(!(chem_class.IsPeptideLinking() || chem_class.IsNucleotideLinking())) {
        throw ost::io::IOException("Could not guess entity type");
      }
    }

    return "polymer";
  }

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

  inline String mon_id_to_olc(char chem_class,
                              const String& mon_id) {

    // hardcoded table according
    // https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entity_poly.pdbx_seq_one_letter_code.html

    if(ost::mol::ChemClass(chem_class).IsPeptideLinking()) {
      switch(mon_id[0]) {
        case 'A': {
          if(mon_id == "ALA") {
            return "A";
          }
          if(mon_id == "ACE") {
            return "(ACE)";
          }
          if(mon_id == "ASP") {
            return "D";
          }
          if(mon_id == "ASN") {
            return "N";
          }
          if(mon_id == "ARG") {
            return "R";
          }
          break;
        }
        case 'C': {
          if(mon_id == "CYS") {
            return "C";
          }
          break;
        }
        case 'G': {
          if(mon_id == "GLU") {
            return "E";
          }
          if(mon_id == "GLY") {
            return "G";
          }
          if(mon_id == "GLN") {
            return "Q";
          }
          break;
        }
        case 'H': {
          if(mon_id == "HIS") {
            return "H";
          }
          break;
        }
        case 'I': {
          if(mon_id == "ILE") {
            return "I";
          }
          break;
        }
        case 'L': {
          if(mon_id == "LEU") {
            return "L";
          }
          if(mon_id == "LYS") {
            return "K";
          }
          break;
        }
        case 'M': {
          if(mon_id == "MET") {
            return "M";
          }
          if(mon_id == "MSE") {
            return "(MSE)";
          }
          break;
        }
        case 'N': {
          if(mon_id == "NH2") {
            return "(NH2)";
          }
          break;
        }
        case 'P': {
          if(mon_id == "PHE") {
            return "F";
          }
          if(mon_id == "PYL") {
            return "O";
          }
          if(mon_id == "PRO") {
            return "P";
          }
          if(mon_id == "PTR") {
            return "(PTR)";
          }
          if(mon_id == "PCA") {
            return "(PCA)";
          }
          break;
        }
        case 'S': {
          if(mon_id == "SER") {
            return "S";
          }
          if(mon_id == "SEC") {
            return "U";
          }
          if(mon_id == "SEP") {
            return "(SEP)";
          }
          break;
        }
        case 'T': {
          if(mon_id == "THR") {
            return "T";
          }
          if(mon_id == "TRP") {
            return "W";
          }
          if(mon_id == "TYR") {
            return "Y";
          }
          if(mon_id == "TPO") {
            return "(TPO)"; 
          }
          break;
        }
        case 'V': {
          if(mon_id == "VAL") {
            return "V";
          }
          break;
        }
      }
    } else if(ost::mol::ChemClass(chem_class).IsNucleotideLinking()) {
      switch(mon_id[0]) {
        case 'A': {
          if(mon_id == "A") {
            return "A";
          }
          break;
        }
        case 'C': {
          if(mon_id == "C") {
            return "C";
          }
          break;
        }
        case 'D': {
          if(mon_id == "DA") {
            return "(DA)";
          }
          if(mon_id == "DC") {
            return "(DC)";
          }
          if(mon_id == "DG") {
            return "(DG)";
          }
          if(mon_id == "DT") {
            return "(DT)";
          }
          break;
        }
        case 'G': {
          if(mon_id == "G") {
            return "G";
          }
          break;
        }
        case 'I': {
          if(mon_id == "I") {
            return "I";
          }
          break;
        } 
        case 'U': {
          if(mon_id == "U") {
            return "U";
          }
          break;
        } 
      }
    } else {
      throw ost::io::IOException("Can only get OLCs for peptides/nucleotides");
    }

    return "(" + mon_id + ")";
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
    String type; // relevant for _entity
    String poly_type; // relevant for _entity_poly
    std::vector<String> asym_ids; // relevant for _struct_asym.id
    std::vector<String> mon_ids; // relevant for _entity_poly_seq.mon_id
    bool is_poly; // in principle type == "polymer"
    String seq; // _entity_poly.pdbx_seq_one_letter_code
    String seq_can; // _entity_poly.pdbx_seq_one_letter_code_can 
  };

  int Setup_entity_(const String& asym_chain_name,
                    const ost::mol::ResidueHandleList& res_list,
                    std::vector<EntityInfo>& entity_infos) {

    String type = GuessEntityType(res_list);
    String poly_type = "";
    bool is_poly = type == "polymer";
    if(is_poly) {
      poly_type = GuessEntityPolyType(res_list);
    }

    std::vector<String> mon_ids;
    if(type == "water") {
      mon_ids.push_back("HOH");
    } else {
      for(auto res : res_list) {
        mon_ids.push_back(res.GetName());
      }
    }

    // check if entity is already there
    for(size_t i = 0; i < entity_infos.size(); ++i) {
      if(entity_infos[i].type == type &&
        entity_infos[i].poly_type == poly_type &&
        entity_infos[i].mon_ids == mon_ids) {
        entity_infos[i].asym_ids.push_back(asym_chain_name);
        return i;
      }
    }

    // need to create new entity
    int entity_idx = entity_infos.size();
    entity_infos.push_back(EntityInfo());
    entity_infos.back().type = type;
    entity_infos.back().poly_type = poly_type;
    entity_infos.back().asym_ids.push_back(asym_chain_name);
    entity_infos.back().mon_ids = mon_ids;
    entity_infos.back().is_poly = is_poly;
    if(is_poly) {
      std::stringstream seq;
      std::stringstream seq_can;
      for(auto res: res_list) {
        // seqres basically follows a hardcoded table from the mmcif reference
        seq << mon_id_to_olc(res.GetChemClass(), res.GetName());
        // canonical seqres maps one letter codes of parent residues 
        // OpenStructure does the same when setting up the entity in the
        // processor. We just trust OpenStructure to do the right thing but set
        // olc to 'X' in case of invalid olc outside [A-Z] (example is '?')
        char olc = res.GetOneLetterCode();
        if(olc < 'A' || olc > 'Z') {
          seq_can << 'X';
        } else {
          seq_can << res.GetOneLetterCode();
        }
      }
      entity_infos.back().seq = seq.str();
      entity_infos.back().seq_can = seq_can.str();
    }

    return entity_idx;
  }

  int Setup_entity_(const String& asym_chain_name,
                    ost::mol::ChainType chain_type,
                    const ost::mol::ResidueHandleList& res_list,
                    std::vector<EntityInfo>& entity_infos) {
    // instead of infering a chain type, we're just checking whether the
    // residues in res_list are allowed to be present in given chain type.
    int entity_idx = 0;
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
    desc.Add("pdbx_seq_one_letter_code");
    desc.Add("pdbx_seq_one_letter_code_can");
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
      ent_data.push_back(ost::io::StarLoopDataItemDO(entity_idx));
      ent_data.push_back(ost::io::StarLoopDataItemDO(entity_info[entity_idx].type));
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
    entity_poly_data.push_back(ost::io::StarLoopDataItemDO("A"));
    entity_poly_data.push_back(ost::io::StarLoopDataItemDO("A"));
    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      if(entity_info[entity_idx].is_poly) {
        entity_poly_data[0] = ost::io::StarLoopDataItemDO(entity_idx);
        entity_poly_data[1] = ost::io::StarLoopDataItemDO(entity_info[entity_idx].poly_type);
        entity_poly_data[2] = ost::io::StarLoopDataItemDO(entity_info[entity_idx].seq);
        entity_poly_data[3] = ost::io::StarLoopDataItemDO(entity_info[entity_idx].seq_can);       
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

  void SetEntity_(const ost::mol::EntityHandle& ent,
                  std::map<String, CompInfo>& comp_infos,
                  std::vector<EntityInfo>& entity_info,
                  ost::io::StarLoop* atom_site,
                  ost::io::StarLoop* pdbx_poly_seq_scheme) {
    ost::mol::ChainHandleList chain_list = ent.GetChainList();
    for(auto ch: chain_list) {

      ost::mol::ResidueHandleList res_list = ch.GetResidueList();

      Setup_chem_comp_(res_list, comp_infos);
      String chain_name = ch.GetName();
      int entity_id = Setup_entity_(chain_name,
                                    ch.GetType(),
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
    }
  }

  void SetEntitymmCIFify_(const ost::mol::EntityHandle& ent,
                          std::map<String, CompInfo>& comp_infos,
                          std::vector<EntityInfo>& entity_info,
                          ost::io::StarLoop* atom_site,
                          ost::io::StarLoop* pdbx_poly_seq_scheme) {

    std::vector<std::vector<ost::mol::ResidueHandle> > L_chains; // L_PEPTIDE_LINKING
    std::vector<std::vector<ost::mol::ResidueHandle> > D_chains; // D_PEPTIDE_LINKING
    std::vector<std::vector<ost::mol::ResidueHandle> > P_chains; // PEPTIDE_LINKING  
    std::vector<std::vector<ost::mol::ResidueHandle> > R_chains; // RNA_LINKING
    std::vector<std::vector<ost::mol::ResidueHandle> > S_chains; // DNA_LINKING
    // all Saccharides go into the same chain
    std::vector<std::vector<ost::mol::ResidueHandle> > Z_chains; // SACCHARIDE 
    std::vector<std::vector<ost::mol::ResidueHandle> > W_chains; // WATER
    std::vector<ost::mol::ResidueHandle> N_chains; // NON_POLYMER (1 res per chain)

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
      bool has_saccharide = false;
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
          has_saccharide = true;
        }

        if(res.GetChemClass() == ost::mol::ChemClass::D_SACCHARIDE) {
          has_saccharide = true;
        }

        if(res.GetChemClass() == ost::mol::ChemClass::SACCHARIDE) {
          has_saccharide = true;
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

      if(has_saccharide) {
        Z_chains.push_back(ost::mol::ResidueHandleList());
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
          Z_chains.back().push_back(res);
        } else if(res.GetChemClass() == ost::mol::ChemClass::D_SACCHARIDE) {
          Z_chains.back().push_back(res);
        } else if(res.GetChemClass() == ost::mol::ChemClass::SACCHARIDE) {
          Z_chains.back().push_back(res);
        } else if(res.GetChemClass() == ost::mol::ChemClass::WATER) {
          W_chains.back().push_back(res);
        } else if(res.GetChemClass() == ost::mol::ChemClass::NON_POLYMER) {
          N_chains.push_back(res);
        } else if(res.GetChemClass() == ost::mol::ChemClass::UNKNOWN) {
          // unknown is just treated as non-poly
          N_chains.push_back(res);
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
      int entity_id = Setup_entity_(chain_name,
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
    }

    // process D_PEPTIDE_LINKING
    for(auto res_list: D_chains) {
      String chain_name = chain_name_gen.Get();
      int entity_id = Setup_entity_(chain_name,
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
    }

    // process PEPTIDE_LINKING
    for(auto res_list: P_chains) {
      String chain_name = chain_name_gen.Get();
      int entity_id = Setup_entity_(chain_name,
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
    }

    // process RNA_LINKING
    for(auto res_list: R_chains) {
      String chain_name = chain_name_gen.Get();
      int entity_id = Setup_entity_(chain_name,
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
    }

    // process DNA_LINKING
    for(auto res_list: S_chains) {
      String chain_name = chain_name_gen.Get();
      int entity_id = Setup_entity_(chain_name,
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
    }

    // process SACHARIDE
    for(auto res_list: Z_chains) {
      String chain_name = chain_name_gen.Get();
      int entity_id = Setup_entity_(chain_name,
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
    }

    // process WATER
    for(auto res_list: W_chains) {
      String chain_name = chain_name_gen.Get();
      int entity_id = Setup_entity_(chain_name,
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
    }

    // process NON_POLYMER
    for(auto res: N_chains) {
      ost::mol::ResidueHandleList res_list;
      res_list.push_back(res);
      String chain_name = chain_name_gen.Get();
      int entity_id = Setup_entity_(chain_name,
                                    res_list,
                                    entity_info);
      Feed_atom_site_(atom_site, chain_name, entity_id, res_list,
                      entity_info[entity_id].is_poly);
      if(entity_info[entity_id].is_poly) {
        Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                  entity_id, res_list);
      }
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

void MMCifWriter::SetEntity(const ost::mol::EntityHandle& ent,
                            bool mmcif_conform) {

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

  std::map<String, CompInfo> comp_infos;
  std::vector<EntityInfo> entity_info;

  // The SetEntity functions fill CompInfo and EntityInfo, i.e. gather info on
  // all unique compounds and entities observed in ent and relate the entities
  // with asym chain names that are directly written to
  // atom_site_/pdbx_poly_seq_scheme_.
  if(mmcif_conform) {
    // chains are assumed to be mmCIF conform - that means water in separate
    // chains, ligands in separate chains etc. Chain types are inferred from
    // chain type property set to the chains in ent. If any of the residues in
    // the respective chains is incompatible with the set chain type, an error
    // is thrown.
    SetEntity_(ent, comp_infos, entity_info,
               atom_site_, pdbx_poly_seq_scheme_);
  } else {
    // rule based splitting of chains into mmCIF conform chains
    SetEntitymmCIFify_(ent, comp_infos, entity_info,
                       atom_site_, pdbx_poly_seq_scheme_);
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
