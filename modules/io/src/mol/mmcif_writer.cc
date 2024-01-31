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

#include <unordered_set>

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

  void CheckValidEntityPolyType(const String& entity_poly_type) {
    std::unordered_set<std::string> s = {"cyclic-pseudo-peptide",
                                         "other",
                                         "peptide nucleic acid",
                                         "polydeoxyribonucleotide",
                                         "polydeoxyribonucleotide/polyribonucleotide hybrid",
                                         "polypeptide(D)",
                                         "polypeptide(L)",
                                         "polyribonucleotide"};
    if(s.find(entity_poly_type) == s.end()) {
      std::stringstream ss;
      ss << "Observed value is no valid entity_poly.type: \"";
      ss << entity_poly_type << "\". Allowed values: ";
      for(auto type: s) {
        ss << type << ", ";
      }
      String err = ss.str();
      throw ost::io::IOException(err.substr(0, err.size() - 2));
    }
  }

  inline String ChemClassToChemCompType(char chem_class) {
    String type = "";
    switch(chem_class) {
      case 'P': {
        type = "PEPTIDE LINKING";
        break;
      }
      case 'D': {
        type = "D-PEPTIDE LINKING";
        break;
      }
      case 'L': {
        type = "L-PEPTIDE LINKING";
        break;
      }
      case 'R': {
        type = "RNA LINKING";
        break;
      }
      case 'S': {
        type = "DNA LINKING";
        break;
      }
      case 'N': {
        type = "NON-POLYMER";
        break;
      }
      case 'X': {
        type = "L-SACCHARIDE";
        break;
      }
      case 'Y': {
        type = "D-SACCHARIDE";
        break;
      }
      case 'Z': {
        type = "SACCHARIDE";
        break;
      }
      case 'W': {
        type = "NON-POLYMER"; // yes, water is a non-polymer
                              // https://www.rcsb.org/ligand/HOH
        break;
      }
      case 'U': {
        type = "OTHER";
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

  inline String MonIDToOLC(const String& mon_id) {

    // hardcoded table according
    // https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entity_poly.pdbx_seq_one_letter_code.html

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
        if(mon_id == "A") {
          return "A";
        }
        break;
      }
      case 'C': {
        if(mon_id == "CYS") {
          return "C";
        }
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
        if(mon_id == "GLU") {
          return "E";
        }
        if(mon_id == "GLY") {
          return "G";
        }
        if(mon_id == "GLN") {
          return "Q";
        }
        if(mon_id == "G") {
          return "G";
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
        if(mon_id == "I") {
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
      case 'U': {
        if(mon_id == "U") {
          return "U";
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

    return "(" + mon_id + ")";
  }

  // template to allow ost::mol::ResidueHandleList and ost::mol::ResidueViewList
  template<class T>
  bool MatchEntity(const T& res_list,
                   const ost::io::MMCifWriterEntity& info) {
    // checks if the residue names in res_list are an exact match
    // with mon_ids in info
    std::vector<String> mon_ids;
    for(auto res : res_list) {
      mon_ids.push_back(res.GetName());
    }
    return mon_ids == info.mon_ids;
  }

  void AddAsym(const String& asym_chain_name,
               ost::io::MMCifWriterEntity& info) {
    // adds asym_chain_name to info under the assumption that mon_ids
    // exactly match => just add a copy of mon_ids to asym_alns
    info.asym_ids.push_back(asym_chain_name);
    info.asym_alns.push_back(info.mon_ids);
  }

  // template to allow ost::mol::ResidueHandleList and ost::mol::ResidueViewList
  template<class T>
  bool MatchEntityResnum(const T& res_list,
                         const ost::io::MMCifWriterEntity& info,
                         Real beyond_frac = 0.05) {
    // Checks if res_list matches SEQRES given in info.mon_ids
    // It may be that res_list defines residues beyond this SEQRES or
    // has residues that are not yet defined, i.e. the respective mon_id is "-".
    // This function returns True if the number of such occurences is below
    // the specified fraction of actual matches (default: 5%).
    int n_beyond = 0;
    int n_matches = 0;
    int num_mon_ids = info.mon_ids.size();
    for(auto res: res_list) {
      int num = res.GetNumber().GetNum();
      char ins_code = res.GetNumber().GetInsCode();
      if(num < 1) {
        std::stringstream ss;
        ss << "Try to construct resnum based alignments. Negative residue ";
        ss << "numbers are not allowed in this case. Got: ";
        ss << num << " in residue " << res;
        ss << ". You may set mmcif_conform flag to False to write something ";
        ss << "but be aware of the consequences...";
        throw ost::io::IOException(ss.str());
      }
      if(ins_code != '\0') {
        std::stringstream ss;
        ss << "Try to construct resnum based alignments. Insertion codes ";
        ss << "are not allowed in this case. Got: ";
        ss << ins_code << " in residue " << res;
        ss << ". You may set mmcif_conform flag to False to write something ";
        ss << "but be aware of the consequences...";
        throw ost::io::IOException(ss.str());
      }
      if (num > num_mon_ids) {
        ++n_beyond;
      } else {
        if(info.mon_ids[num-1] == "-") {
          ++n_beyond; // we're basically filling an unknown gap...
        } else if(info.mon_ids[num-1] == res.GetName()) {
          ++n_matches;
        } else {
          return false;
        }
      }
    }
    if(n_matches == 0) {
      return false;
    } else {
      return n_beyond / n_matches <= beyond_frac;
    }
  }

  // template to allow ost::mol::ResidueHandleList and ost::mol::ResidueViewList
  template<class T>
  void AddAsymResnum(const String& asym_chain_name,
                     const T& res_list,
                     ost::io::MMCifWriterEntity& info) {

    if(!info.is_poly) {
      // no need for SEQRES alignment vodoo
      AddAsym(asym_chain_name, info);
      return;
    }

    int max_resnum = info.mon_ids.size();
    std::vector<String> mon_ids;
    std::vector<int> resnums;

    for(auto res: res_list) {
      int num = res.GetNumber().GetNum();
      // assumes that MatchEntityResnum has already been run to check for
      // resnum < 1
      mon_ids.push_back(res.GetName());
      resnums.push_back(num);
      max_resnum = std::max(max_resnum, num);
    }

    std::vector<String> aln_mon_ids(max_resnum, "-");
    for(size_t i = 0; i < mon_ids.size(); ++i) {
      aln_mon_ids[resnums[i]-1] = mon_ids[i];
    }

    if(max_resnum > static_cast<int>(info.mon_ids.size())) {
      // This chain covers more residues towards C-terminus than any chain that
      // is associated with this entity - expand to enforce equal size
      int N = max_resnum - info.mon_ids.size();
      info.mon_ids.insert(info.mon_ids.end(), N, "-");
      info.seq_olcs.insert(info.seq_olcs.end(), N, "-");
      info.seq_can_olcs.insert(info.seq_can_olcs.end(), N, "-");
      for(size_t asym_idx = 0; asym_idx < info.asym_alns.size(); ++asym_idx) {
        info.asym_alns[asym_idx].insert(info.asym_alns[asym_idx].end(), N, "-");
      }
    }

    // Fill SEQRES infos newly covered by this asym chain
    for(size_t i = 0; i < resnums.size(); ++i) {
      if(info.mon_ids[resnums[i]-1] == "-") {
        info.mon_ids[resnums[i]-1] = mon_ids[i];
        info.seq_olcs[resnums[i]-1] = MonIDToOLC(mon_ids[i]);
        char olc = res_list[i].GetOneLetterCode();
        if(olc < 'A' || olc > 'Z') {
          info.seq_can_olcs[resnums[i]-1] = "X";
        } else {
          info.seq_can_olcs[resnums[i]-1] = String(1, olc);
        }
      }
    }
    
    // finalize
    info.asym_ids.push_back(asym_chain_name);
    info.asym_alns.push_back(aln_mon_ids);
  }

  // template to allow ost::mol::ResidueHandleList and ost::mol::ResidueViewList
  template<class T>
  int SetupEntity(const String& asym_chain_name,
                  const String& type,
                  const String& poly_type,
                  const String& branch_type,
                  const T& res_list,
                  bool resnum_alignment,
                  std::vector<ost::io::MMCifWriterEntity>& entity_infos) {

    bool is_poly = type == "polymer";

    if(!is_poly && res_list.size() != 1 && type != "water" && type != "branched") {
      std::stringstream ss;
      ss << "Cannot setup entity with " << res_list.size() << " residues ";
      ss << "but is of type: " << type;
      throw ost::io::IOException(ss.str());
    }

    // check if entity is already there
    for(size_t i = 0; i < entity_infos.size(); ++i) {
      if(entity_infos[i].type == "water" && type == "water") {
        AddAsym(asym_chain_name, entity_infos[i]);
        return i;
      }
      if(entity_infos[i].type == type &&
         entity_infos[i].poly_type == poly_type &&
         entity_infos[i].branch_type == branch_type) {
        if(is_poly && resnum_alignment) {
          if(MatchEntityResnum(res_list, entity_infos[i])) {
            AddAsymResnum(asym_chain_name, res_list, entity_infos[i]);
            return i;
          }   
        } else {
           if(MatchEntity(res_list, entity_infos[i])) {
             AddAsym(asym_chain_name, entity_infos[i]);
             return i;
           }
        }
      }
    }

    // need to create new entity
    std::vector<String> mon_ids;
    std::vector<String> seq;
    std::vector<String> seq_can;

    if(is_poly && resnum_alignment) {
      int max_resnum = res_list.size();
      std::vector<String> res_mon_ids;
      std::vector<int> resnums;
      for(auto res: res_list) {
        int num = res.GetNumber().GetNum();
        char ins_code = res.GetNumber().GetInsCode();
        if(num < 1) {
          std::stringstream ss;
          ss << "Try to construct mmCIF entity from residues using resnum ";
          ss << "based alignments. Negative residue numbers are not allowed ";
          ss << "in this case. Got: " << num << " in residue " << res;
          ss << ". You may set mmcif_conform flag to False to write something ";
          ss << "but be aware of the consequences...";
          throw ost::io::IOException(ss.str());
        }
        if(ins_code != '\0') {
          std::stringstream ss;
          ss << "Try to construct mmCIF entity from residues using resnum ";
          ss << "based alignments. Insertion codes are not allowed ";
          ss << "in this case. Got: " << ins_code << " in residue " << res;
          ss << ". You may set mmcif_conform flag to False to write something ";
          ss << "but be aware of the consequences...";
          throw ost::io::IOException(ss.str());
        }
        res_mon_ids.push_back(res.GetName());
        resnums.push_back(num);
        max_resnum = std::max(max_resnum, num);
      }
      mon_ids.assign(max_resnum, "-");
      seq.assign(max_resnum, "-");
      seq_can.assign(max_resnum, "-");
      for(size_t i = 0; i < res_mon_ids.size(); ++i) {
        mon_ids[resnums[i]-1] = res_mon_ids[i];
        seq[resnums[i]-1] = MonIDToOLC(mon_ids[resnums[i]-1]);
        char olc = res_list[i].GetOneLetterCode();
        if(olc < 'A' || olc > 'Z') {
          seq_can[resnums[i]-1] = "X";
        } else {
          seq_can[resnums[i]-1] = String(1, olc);
        }
      }
    } else {
      if(type == "water") {
        mon_ids.push_back("HOH");
        seq.push_back("(HOH)");
        seq_can.push_back("?");
      } else {
        for(auto res: res_list) {
          mon_ids.push_back(res.GetName());
          seq.push_back(MonIDToOLC(res.GetName()));
          char olc = res.GetOneLetterCode();
          if(olc < 'A' || olc > 'Z') {
            seq_can.push_back("X");
          } else {
            seq_can.push_back(String(1, olc));
          }
        }
      }
    }

    int entity_idx = entity_infos.size();
    entity_infos.push_back(ost::io::MMCifWriterEntity());
    entity_infos.back().type = type;
    entity_infos.back().poly_type = poly_type;
    entity_infos.back().branch_type = branch_type;
    entity_infos.back().mon_ids = mon_ids;
    entity_infos.back().seq_olcs = seq;
    entity_infos.back().seq_can_olcs = seq_can;
    entity_infos.back().is_poly = is_poly;

    if(is_poly && resnum_alignment) {
      AddAsymResnum(asym_chain_name, res_list, entity_infos.back());
    } else {
      AddAsym(asym_chain_name, entity_infos.back());
    }

    return entity_idx;
  }

  // template to allow ost::mol::ResidueHandleList and ost::mol::ResidueViewList
  template<class T>
  int SetupEntity(const String& asym_chain_name,
                  ost::mol::ChainType chain_type,
                  const T& res_list,
                  bool resnum_alignment, 
                  std::vector<ost::io::MMCifWriterEntity>& entity_infos) {
    // use chain_type info attached to chain to determine
    // _entity.type and _entity_poly.type
    String type = ost::mol::EntityTypeFromChainType(chain_type);
    bool is_poly = type == "polymer";
    String poly_type = "";
    if(is_poly) {
      poly_type = ost::mol::EntityPolyTypeFromChainType(chain_type);
    }
    bool is_branched = type == "branched";
    String branch_type = "";
    if(is_branched) {
      branch_type = ost::mol::BranchedTypeFromChainType(chain_type);
    }
    return SetupEntity(asym_chain_name, type, poly_type, branch_type, res_list,
                       resnum_alignment, entity_infos);
  }

  ost::io::StarWriterLoopPtr Setup_atom_type_ptr() {
    ost::io::StarWriterLoopDesc desc("_atom_type");
    desc.Add("symbol");
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;
  }

  ost::io::StarWriterLoopPtr Setup_atom_site_ptr() {
    ost::io::StarWriterLoopDesc desc("_atom_site");
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
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;
  }

  ost::io::StarWriterLoopPtr Setup_pdbx_poly_seq_scheme_ptr() {
    ost::io::StarWriterLoopDesc desc("_pdbx_poly_seq_scheme");
    desc.Add("asym_id");
    desc.Add("entity_id");
    desc.Add("mon_id");
    desc.Add("seq_id");
    desc.Add("pdb_strand_id");
    desc.Add("pdb_seq_num");
    desc.Add("pdb_ins_code");
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;
  }

  ost::io::StarWriterLoopPtr Setup_entity_ptr() {
    ost::io::StarWriterLoopDesc desc("_entity");
    desc.Add("id");
    desc.Add("type");
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;
  }

  ost::io::StarWriterLoopPtr Setup_struct_asym_ptr() {
    ost::io::StarWriterLoopDesc desc("_struct_asym");
    desc.Add("id");
    desc.Add("entity_id");
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;    
  }

  ost::io::StarWriterLoopPtr Setup_entity_poly_ptr() {
    ost::io::StarWriterLoopDesc desc("_entity_poly");
    desc.Add("entity_id");
    desc.Add("type");
    desc.Add("pdbx_seq_one_letter_code");
    desc.Add("pdbx_seq_one_letter_code_can");
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;    
  }

  ost::io::StarWriterLoopPtr Setup_entity_poly_seq_ptr() {
    ost::io::StarWriterLoopDesc desc("_entity_poly_seq");
    desc.Add("entity_id");
    desc.Add("mon_id");
    desc.Add("num");
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;    
  }

  ost::io::StarWriterLoopPtr Setup_chem_comp_ptr() {
    ost::io::StarWriterLoopDesc desc("_chem_comp");
    desc.Add("id");
    desc.Add("type");
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;    
  }

  ost::io::StarWriterLoopPtr Setup_pdbx_entity_branch_ptr() {
    ost::io::StarWriterLoopDesc desc("_pdbx_entity_branch");
    desc.Add("entity_id");
    desc.Add("type");
    ost::io::StarWriterLoopPtr sl(new ost::io::StarWriterLoop(desc));
    return sl;    
  }

  void Feed_atom_type(ost::io::StarWriterLoopPtr atom_type_ptr,
                      ost::io::StarWriterLoopPtr atom_site_ptr) {
    // we're just extracting every type_symbol that we observed
    // in atom_site (this is a bit of circular stupidity...)
    std::set<String> symbols;
    int desc_size = atom_site_ptr->GetDesc().GetSize();
    int type_symbol_idx = atom_site_ptr->GetDesc().GetIndex("type_symbol");
    int N = atom_site_ptr->GetN();
    const std::vector<ost::io::StarWriterValue>& data = atom_site_ptr->GetData();
    for(int i = 0; i < N; ++i) {
      symbols.insert(data[i*desc_size + type_symbol_idx].GetValue());
    }
    std::vector<ost::io::StarWriterValue> atom_type_data(1);
    for(auto symbol: symbols) {
      atom_type_data[0] = ost::io::StarWriterValue::FromString(symbol);
      atom_type_ptr->AddData(atom_type_data);
    }
  }

  // template to allow ost::mol::ResidueHandleList and ost::mol::ResidueViewList
  template<class T>
  void Feed_pdbx_poly_seq_scheme(ost::io::StarWriterLoopPtr pdbx_poly_seq_scheme_ptr,
                                 const String& label_asym_id,
                                 int label_entity_id,
                                 ost::io::MMCifWriterEntity& entity_info,
                                 const T& res_list) {

    std::vector<ost::io::StarWriterValue> data(7);
    // processing chain by chain, label_asym_id and label_entity_id are constant
    data[0] = ost::io::StarWriterValue::FromString(label_asym_id);
    data[1] = ost::io::StarWriterValue::FromInt(label_entity_id);

    int asym_idx = entity_info.GetAsymIdx(label_asym_id);
    const std::vector<String>& aln = entity_info.asym_alns[asym_idx];
    int label_seq_id = 0; // 0-based index

    for(auto res: res_list) {
      String res_name = res.GetName();
      while(aln[label_seq_id] == "-") {
        ++label_seq_id;
      }

      data[2] = ost::io::StarWriterValue::FromString(res_name);
      data[3] = ost::io::StarWriterValue::FromInt(label_seq_id + 1);

      // the remaining data items honor String properties if set:
      // pdb_auth_chain_name, pdb_auth_resnum and pdb_auth_ins_code

      if(res.GetChain().HasProp("pdb_auth_chain_name")) {
        data[4] = 
        ost::io::StarWriterValue::FromString(res.GetChain().GetStringProp("pdb_auth_chain_name"));
      } else {
        data[4] = ost::io::StarWriterValue::FromString(res.GetChain().GetName());  
      }

      if(res.HasProp("pdb_auth_resnum")) {
        // this feels so wrong that this is stored as a string property...
        data[5] = ost::io::StarWriterValue::FromString(res.GetStringProp("pdb_auth_resnum"));
      } else {
        data[5] = ost::io::StarWriterValue::FromInt(res.GetNumber().GetNum());
      }

      if(res.HasProp("pdb_auth_ins_code")) {
        data[6] = ost::io::StarWriterValue::FromString(res.GetStringProp("pdb_auth_ins_code"));
      } else {
        char ins_code = res.GetNumber().GetInsCode();      
        if(ins_code == '\0') {
          data[6] = ost::io::StarWriterValue::FromString("");
        } else {
          data[6] = ost::io::StarWriterValue::FromString(String(1, ins_code));
        }      
      }
      pdbx_poly_seq_scheme_ptr->AddData(data);
      label_seq_id += 1;
    }
  }

  // template to allow ost::mol::ResidueHandleList and ost::mol::ResidueViewList
  template<class T>
  void Feed_atom_site(ost::io::StarWriterLoopPtr atom_site_ptr,
                      const String& label_asym_id,
                      int label_entity_id,
                      const ost::io::MMCifWriterEntity& entity_info,
                      const T& res_list) {

    int asym_idx = entity_info.GetAsymIdx(label_asym_id);
    const std::vector<String>& aln = entity_info.asym_alns[asym_idx];
    int label_seq_id = 0; // 0-based index
    std::vector<ost::io::StarWriterValue> at_data(17);

    for(auto res: res_list) {
      String comp_id = res.GetName();

      auto at_list = res.GetAtomList();
      String auth_asym_id = res.GetChain().GetName();
      if(res.HasProp("pdb_auth_chain_name")) {
        auth_asym_id = res.GetStringProp("pdb_auth_chain_name");
      }
      
      String auth_seq_id = std::to_string(res.GetNumber().GetNum());
      if(res.HasProp("pdb_auth_resnum")) {
        auth_seq_id = res.GetStringProp("pdb_auth_resnum");
      }

      char c_ins_code = res.GetNumber().GetInsCode();
      String ins_code = c_ins_code == '\0' ? "" : String(1, c_ins_code);
      if(res.HasProp("pdb_auth_ins_code")) {
        ins_code = res.GetStringProp("pdb_auth_ins_code");
      }

      if(entity_info.is_poly) {
        while(aln[label_seq_id] == "-") {
          ++label_seq_id;
        }
      }

      for(auto at: at_list) {
        // group_PDB
        if(at.IsHetAtom()) {
          at_data[0] = ost::io::StarWriterValue::FromString("HETATM");
        } else {
          at_data[0] = ost::io::StarWriterValue::FromString("ATOM");
        }
        // type_symbol
        at_data[1] = ost::io::StarWriterValue::FromString(at.GetElement());
        // label_atom_id
        at_data[2] = ost::io::StarWriterValue::FromString(at.GetName());
        // label_comp_id
        at_data[3] = ost::io::StarWriterValue::FromString(comp_id);
        // label_asym_id
        at_data[4] = ost::io::StarWriterValue::FromString(label_asym_id);
        // label_entity_id
        at_data[5] = ost::io::StarWriterValue::FromInt(label_entity_id);
        // label_seq_id
        if(entity_info.is_poly) {
          at_data[6] = ost::io::StarWriterValue::FromInt(label_seq_id+1);
        } else {
          at_data[6] = ost::io::StarWriterValue::FromString(".");
        }
        // label_alt_id
        at_data[7] = ost::io::StarWriterValue::FromString(".");
        // Cartn_x
        at_data[8] = ost::io::StarWriterValue::FromFloat(at.GetPos().GetX(), 3);
        // Cartn_y
        at_data[9] = ost::io::StarWriterValue::FromFloat(at.GetPos().GetY(), 3);
        // Cartn_z
        at_data[10] = ost::io::StarWriterValue::FromFloat(at.GetPos().GetZ(), 3);
        // occupancy
        at_data[11] = ost::io::StarWriterValue::FromFloat(at.GetOccupancy(), 2);
        // B_iso_or_equiv
        at_data[12] = ost::io::StarWriterValue::FromFloat(at.GetBFactor(), 2);
        // auth_seq_id
        at_data[13] = ost::io::StarWriterValue::FromString(auth_seq_id);
        // auth_asym_id
        at_data[14] = ost::io::StarWriterValue::FromString(auth_asym_id);
        // id
        at_data[15] = ost::io::StarWriterValue::FromInt(atom_site_ptr->GetN());
        // pdbx_PDB_ins_code
        at_data[16] = ost::io::StarWriterValue::FromString(ins_code);
        atom_site_ptr->AddData(at_data);
      }
      ++label_seq_id;
    }
  }

  void Feed_entity(ost::io::StarWriterLoopPtr entity_ptr,
                   const std::vector<ost::io::MMCifWriterEntity>& entity_info) {
    std::vector<ost::io::StarWriterValue> ent_data(2);
    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      ent_data[0] = ost::io::StarWriterValue::FromInt(entity_idx);
      ent_data[1] = ost::io::StarWriterValue::FromString(entity_info[entity_idx].type);
      entity_ptr->AddData(ent_data);
    }
  }

  void Feed_struct_asym(ost::io::StarWriterLoopPtr struct_asym_ptr,
                        const std::vector<ost::io::MMCifWriterEntity>& entity_info) {
    std::vector<ost::io::StarWriterValue> asym_data(2);
    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      for(auto asym_id : entity_info[entity_idx].asym_ids) {
        asym_data[0] = ost::io::StarWriterValue::FromString(asym_id);
        asym_data[1] = ost::io::StarWriterValue::FromInt(entity_idx);
        struct_asym_ptr->AddData(asym_data);
      }
    }
  }

  void Feed_entity_poly_seq(ost::io::StarWriterLoopPtr entity_poly_seq_ptr,
                            const std::vector<ost::io::MMCifWriterEntity>& entity_info) {
    std::vector<ost::io::StarWriterValue> entity_poly_seq_data(3);
    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      if(entity_info[entity_idx].is_poly) {
        const std::vector<String>& mon_ids = entity_info[entity_idx].mon_ids;
        for(size_t mon_idx = 0; mon_idx < mon_ids.size(); ++mon_idx) {
          entity_poly_seq_data[0] = ost::io::StarWriterValue::FromInt(entity_idx);
          entity_poly_seq_data[1] = ost::io::StarWriterValue::FromString(mon_ids[mon_idx]);
          entity_poly_seq_data[2] = ost::io::StarWriterValue::FromInt(mon_idx+1);
          entity_poly_seq_ptr->AddData(entity_poly_seq_data);
        }
      }
    }
  }

  void Feed_entity_poly(ost::io::StarWriterLoopPtr entity_poly_ptr,
                        const std::vector<ost::io::MMCifWriterEntity>& entity_info) {
    std::vector<ost::io::StarWriterValue> entity_poly_data(4);
    for(size_t entity_idx = 0; entity_idx < entity_info.size(); ++entity_idx) {
      if(entity_info[entity_idx].is_poly) {
        entity_poly_data[0] = ost::io::StarWriterValue::FromInt(entity_idx);
        entity_poly_data[1] = ost::io::StarWriterValue::FromString(entity_info[entity_idx].poly_type);
        std::stringstream seq;
        std::stringstream seq_can;
        for(size_t idx = 0; idx < entity_info[entity_idx].mon_ids.size(); ++idx) {
          if(entity_info[entity_idx].seq_olcs[idx] == "-") {
            // X IS PROBABLY NOT THE RIGHT THING FOR
            // pdbx_seq_one_letter_code
            seq << "X";
          } else {
            seq << entity_info[entity_idx].seq_olcs[idx];
          }
          if(entity_info[entity_idx].seq_can_olcs[idx] == "-") {
            seq_can << "X";
          } else {
            seq_can << entity_info[entity_idx].seq_can_olcs[idx];
          }
        }
        entity_poly_data[2] = ost::io::StarWriterValue::FromString(seq.str());
        entity_poly_data[3] = ost::io::StarWriterValue::FromString(seq_can.str());
        entity_poly_ptr->AddData(entity_poly_data);
      }
    }
  }

  void Feed_chem_comp(ost::io::StarWriterLoopPtr chem_comp_ptr,
                      const std::vector<ost::io::MMCifWriterEntity>& entity_infos,
                      ost::conop::CompoundLibPtr compound_lib) {
    std::set<String> unique_compounds;
    for(auto ent: entity_infos) {
      unique_compounds.insert(ent.mon_ids.begin(), ent.mon_ids.end());
    }
    std::vector<ost::io::StarWriterValue> comp_data(2);
    for(auto mon_id: unique_compounds) {
      comp_data[0] = ost::io::StarWriterValue::FromString(mon_id);
      ost::conop::CompoundPtr comp = compound_lib->FindCompound(mon_id,
                                                                ost::conop::Compound::PDB);
      if(comp) {
        String type = ChemClassToChemCompType(comp->GetChemClass());
        comp_data[1] = ost::io::StarWriterValue::FromString(type);
      } else {
        String type = ChemClassToChemCompType(ost::mol::ChemClass::UNKNOWN);
        comp_data[1] = ost::io::StarWriterValue::FromString(type);
      }
      chem_comp_ptr->AddData(comp_data);
    }
  }

  void Feed_pdbx_entity_branch(ost::io::StarWriterLoopPtr pdbx_entity_branch_ptr,
                               const std::vector<ost::io::MMCifWriterEntity>& entity_infos) {
    std::vector<ost::io::StarWriterValue> branch_data(2);
    for(size_t i = 0; i < entity_infos.size(); ++i) {
      if(entity_infos[i].type == "branched") {
        branch_data[0] = ost::io::StarWriterValue::FromInt(i);
        branch_data[1] = ost::io::StarWriterValue::FromString(entity_infos[i].branch_type);
        pdbx_entity_branch_ptr->AddData(branch_data);
      }
    }
  }

  // template to allow ost::mol::ResidueHandleList and ost::mol::ResidueViewList
  template<class T>
  void ProcessEntmmCIFify(const std::vector<T>& res_lists,
                          ost::conop::CompoundLibPtr compound_lib,
                          std::vector<ost::io::MMCifWriterEntity>& entity_info,
                          ost::io::StarWriterLoopPtr atom_site,
                          ost::io::StarWriterLoopPtr pdbx_poly_seq_scheme) {

    ChainNameGenerator chain_name_gen;

    std::set<String> unique_compounds;
    for(auto res_list: res_lists) {
      for(auto res: res_list) {
        unique_compounds.insert(res.GetName());
      }
    }
    std::map<String, ost::mol::ChemClass> chem_class_mapping;
    for(auto mon_id: unique_compounds) {
      ost::conop::CompoundPtr comp = compound_lib->FindCompound(mon_id,
                                                                ost::conop::Compound::PDB);
      if(comp) {
        chem_class_mapping[mon_id] = comp->GetChemClass();
      } else {
        chem_class_mapping[mon_id] = ost::mol::ChemClass(ost::mol::ChemClass::UNKNOWN);
      }
    }

    for(auto res_list: res_lists) {

      T L_chain;
      T D_chain;
      T P_chain;
      T R_chain;
      T S_chain;
      T Z_chain;
      T W_chain;

      std::vector<ost::mol::ChemClass> chem_classes;
      chem_classes.reserve(res_list.size());
      for(auto res: res_list) {
        chem_classes.push_back(chem_class_mapping[res.GetName()]);
      }

      // first scan only concerning peptides...
      // Avoid mix of both in same chain: L-peptide linking, D-peptide linking
      // But we still want to know that in advance as we assign non chiral
      // peptides to either of those 
      bool has_l_peptide_linking = false;
      bool has_d_peptide_linking = false;
      for(auto chem_class: chem_classes) {
        if(chem_class == ost::mol::ChemClass::D_PEPTIDE_LINKING) {
          if(has_l_peptide_linking) {
            throw ost::io::IOException("Cannot write mmCIF when same chain "
                                       "contains D- and L-peptides");
          }
          has_d_peptide_linking = true;
        }
        if(chem_class == ost::mol::ChemClass::L_PEPTIDE_LINKING) {
          if(has_d_peptide_linking) {
            throw ost::io::IOException("Cannot write mmCIF when same chain "
                                       "contains D- and L-peptides");
          }
          has_l_peptide_linking = true;
        }        
      }

      for(size_t i = 0; i < res_list.size(); ++i) {
        if(chem_classes[i].IsPeptideLinking()) {
          if(has_l_peptide_linking) {
            L_chain.push_back(res_list[i]);
          } else if(has_d_peptide_linking) {
            D_chain.push_back(res_list[i]);
          } else {
            P_chain.push_back(res_list[i]);
          }
        } else if(chem_classes[i] == ost::mol::ChemClass::RNA_LINKING) {
          R_chain.push_back(res_list[i]);
        } else if(chem_classes[i] == ost::mol::ChemClass::DNA_LINKING) {
          S_chain.push_back(res_list[i]);
        } else if(chem_classes[i] == ost::mol::ChemClass::L_SACCHARIDE) {
          Z_chain.push_back(res_list[i]);
        } else if(chem_classes[i] == ost::mol::ChemClass::D_SACCHARIDE) {
          Z_chain.push_back(res_list[i]);
        } else if(chem_classes[i] == ost::mol::ChemClass::SACCHARIDE) {
          Z_chain.push_back(res_list[i]);
        } else if(chem_classes[i] == ost::mol::ChemClass::WATER) {
          W_chain.push_back(res_list[i]);
        } else if(chem_classes[i] == ost::mol::ChemClass::NON_POLYMER ||
                  chem_classes[i] == ost::mol::ChemClass::UNKNOWN) {
          // already process non-poly and unknown
          String type = "non-polymer";
          String poly_type = "";
          String branch_type = "";
          T tmp;
          tmp.push_back(res_list[i]);
          String chain_name = chain_name_gen.Get();
          int entity_id = SetupEntity(chain_name,
                                      type,
                                      poly_type,
                                      branch_type,
                                      tmp,
                                      false,
                                      entity_info);
          Feed_atom_site(atom_site, chain_name, entity_id, entity_info[entity_id],
                         tmp);
        } else {
          // this should not happen...
          std::stringstream ss;
          ss << "Unsupported chem class (" << res_list[i].GetChemClass();
          ss << ") for residue "<< res_list[i];
          throw ost::io::IOException(ss.str());
        }
      }
      // process poly chains
      T* poly_chains[5] = {&L_chain, &D_chain, &P_chain, &R_chain, &S_chain};
      String poly_types[5] = {"polypeptide(L)", "polypeptide(D)",
                              "polypeptide(L)", "polyribonucleotide",
                              "polydeoxyribonucleotide"};

      for(int i = 0; i < 5; ++i) {
        if(!poly_chains[i]->empty()) {
          String type = "polymer";
          String poly_type = poly_types[i];
          if(poly_chains[i]->size() <= 2) {
            // must have length of at least 3 to be polymer
            type = "non-polymer";
            poly_type = "";
          }
          String branch_type = "";
          String chain_name = chain_name_gen.Get();
          int entity_id = SetupEntity(chain_name,
                                      type,
                                      poly_type,
                                      branch_type,
                                      *poly_chains[i],
                                      false,
                                      entity_info);
          Feed_atom_site(atom_site, chain_name, entity_id, entity_info[entity_id],
                         *poly_chains[i]);
          // still check whether we're dealing with poly here, we could have a
          // lonely amino acid that ends up as non-poly and doesn't need
          // pdbx_poly_seq_scheme
          if(entity_info[entity_id].is_poly) {
            Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                      entity_id, entity_info[entity_id], *poly_chains[i]);
          }
        }
      }

      // process water chain
      if(!W_chain.empty()) {
        String type = "water";
        String poly_type = "";
        String branch_type = "";
        String chain_name = chain_name_gen.Get();
        int entity_id = SetupEntity(chain_name,
                                    type,
                                    poly_type,
                                    branch_type,
                                    W_chain,
                                    false,
                                    entity_info);
        Feed_atom_site(atom_site, chain_name, entity_id, entity_info[entity_id],
                       W_chain);
      }
      // process sugar chain
      if(!Z_chain.empty()) {
        String type = "branched";
        String poly_type = "";
        String branch_type = "oligosaccharide";
        if(Z_chain.size() == 1) {
          // not really branched if the poor little Zueckerli is alone...
          type = "non-polymer"; 
          branch_type = "";
        }

        String chain_name = chain_name_gen.Get();
        int entity_id = SetupEntity(chain_name,
                                    type,
                                    poly_type,
                                    branch_type,
                                    Z_chain,
                                    false,
                                    entity_info);
        Feed_atom_site(atom_site, chain_name, entity_id, entity_info[entity_id],
                       Z_chain);
      }
    }
  }

  void ProcessEntmmCIFify(const ost::mol::EntityHandle& ent,
                          ost::conop::CompoundLibPtr compound_lib,
                          std::vector<ost::io::MMCifWriterEntity>& entity_info,
                          ost::io::StarWriterLoopPtr atom_site,
                          ost::io::StarWriterLoopPtr pdbx_poly_seq_scheme) {
    std::vector<ost::mol::ResidueHandleList> res_lists;
    ost::mol::ChainHandleList chain_list = ent.GetChainList();
    for(auto ch: chain_list) {
      res_lists.push_back(ch.GetResidueList());
    }
    ProcessEntmmCIFify(res_lists, compound_lib, entity_info,
                       atom_site, pdbx_poly_seq_scheme);
  }

  void ProcessEntmmCIFify(const ost::mol::EntityView& ent,
                          ost::conop::CompoundLibPtr compound_lib,
                          std::vector<ost::io::MMCifWriterEntity>& entity_info,
                          ost::io::StarWriterLoopPtr atom_site,
                          ost::io::StarWriterLoopPtr pdbx_poly_seq_scheme) {
    std::vector<ost::mol::ResidueViewList> res_lists;
    ost::mol::ChainViewList chain_list = ent.GetChainList();
    for(auto ch: chain_list) {
      res_lists.push_back(ch.GetResidueList());
    }
    ProcessEntmmCIFify(res_lists, compound_lib, entity_info,
                       atom_site, pdbx_poly_seq_scheme);
  }

  // template to allow ost::mol::EntityHandle and ost::mol::EntityView
  template<class T>
  void ProcessEnt(const T& ent,
                  ost::conop::CompoundLibPtr compound_lib,
                  bool mmcif_conform,
                  std::vector<ost::io::MMCifWriterEntity>& entity_info,
                  ost::io::StarWriterLoopPtr atom_site,
                  ost::io::StarWriterLoopPtr pdbx_poly_seq_scheme) {

    if(mmcif_conform) {
      auto chain_list = ent.GetChainList();
      for(auto ch: chain_list) {
        auto res_list = ch.GetResidueList();
        String chain_name = ch.GetName();
        int entity_id = SetupEntity(chain_name,
                                    ch.GetType(),
                                    res_list,
                                    true,
                                    entity_info);
        Feed_atom_site(atom_site, chain_name, entity_id, entity_info[entity_id],
                       res_list);
        if(entity_info[entity_id].is_poly) {
          Feed_pdbx_poly_seq_scheme(pdbx_poly_seq_scheme, chain_name,
                                    entity_id, entity_info[entity_id], res_list);
        }
      }
    } else {
      // delegate to more complex ProcessEntmmCIFify
      ProcessEntmmCIFify(ent, compound_lib, entity_info, atom_site,
                         pdbx_poly_seq_scheme);
    }
  }

  void ProcessUnknowns(std::vector<ost::io::MMCifWriterEntity>& entity_infos) {

    for(size_t entity_idx = 0; entity_idx < entity_infos.size(); ++entity_idx) {
      if(entity_infos[entity_idx].is_poly) {
        // scan for "-" in mon_ids
        for(size_t mon_id_idx = 0;
            mon_id_idx < entity_infos[entity_idx].mon_ids.size(); ++mon_id_idx) {

          if(entity_infos[entity_idx].mon_ids[mon_id_idx] == "-") {

            if(entity_infos[entity_idx].poly_type == "polypeptide(D)" ||
               entity_infos[entity_idx].poly_type == "polypeptide(L)" ||
               entity_infos[entity_idx].poly_type == "cyclic-pseudo-peptide" ||
               entity_infos[entity_idx].poly_type == "peptide nucleic acid") {
              entity_infos[entity_idx].mon_ids[mon_id_idx] = "UNK"; 
              entity_infos[entity_idx].seq_olcs[mon_id_idx] = "(UNK)"; 
              entity_infos[entity_idx].seq_can_olcs[mon_id_idx] = "X";
            }

            else if(entity_infos[entity_idx].poly_type == "polydeoxyribonucleotide") {
              entity_infos[entity_idx].mon_ids[mon_id_idx] = "DN"; 
              entity_infos[entity_idx].seq_olcs[mon_id_idx] = "(DN)"; 
              entity_infos[entity_idx].seq_can_olcs[mon_id_idx] = "N";
            }

            else if(entity_infos[entity_idx].poly_type == "polyribonucleotide" ||
               entity_infos[entity_idx].poly_type == "polydeoxyribonucleotide/polyribonucleotide hybrid") {
              entity_infos[entity_idx].mon_ids[mon_id_idx] = "N"; 
              entity_infos[entity_idx].seq_olcs[mon_id_idx] = "N"; 
              entity_infos[entity_idx].seq_can_olcs[mon_id_idx] = "N";
            }

            else {
              std::stringstream ss;
              ss << "Gaps are not supported for polymer chains of type ";
              ss << entity_infos[entity_idx].poly_type;
              throw ost::io::IOException(ss.str());
            }
          } 
        }
      }
    }
  }

} // ns

namespace ost { namespace io {

MMCifWriterEntity MMCifWriterEntity::FromPolymer(const String& entity_poly_type,
                                                 const std::vector<String>& mon_ids,
                                                 conop::CompoundLibPtr compound_lib) {
  CheckValidEntityPolyType(entity_poly_type);
  MMCifWriterEntity ent;
  ent.type = "polymer";
  ent.is_poly = true;
  ent.poly_type = entity_poly_type;
  ent.branch_type = "";
  ent.mon_ids = mon_ids;
  for(auto mon_id: mon_ids) {
    ent.seq_olcs.push_back(MonIDToOLC(mon_id));
    if(ent.seq_olcs.back().size() == 1) {
      ent.seq_can_olcs.push_back(ent.seq_olcs.back());
    } else {
      ost::conop::CompoundPtr compound = 
      compound_lib->FindCompound(mon_id, ost::conop::Compound::PDB);
      char olc = compound->GetOneLetterCode();
      if(olc < 'A' || olc > 'Z') {
        ent.seq_can_olcs.push_back("(" + mon_id + ")");  
      } else {
        ent.seq_can_olcs.push_back(String(1, olc));
      }
    }
  }
  return ent;
}

int MMCifWriterEntity::GetAsymIdx(const String& asym_id) const {
  for(size_t i = 0; i < asym_ids.size(); ++i) {
    if(asym_ids[i] == asym_id) {
      return i;
    }
  }
  std::stringstream ss;
  ss << "Tried to find asym id " << asym_id << "in entity that has only ";
  ss << "the following asym ids: ";
  for(auto i: asym_ids) {
    ss << i << ", ";
  }
  String err = ss.str();
  err = err.substr(0, err.size()-2); // remove last ", "
  throw ost::io::IOException(err);  
}

void MMCifWriter::SetStructure(const ost::mol::EntityHandle& ent,
                               conop::CompoundLibPtr compound_lib,
                               bool mmcif_conform,
                               const std::vector<MMCifWriterEntity>& entity_info) {
  this->Setup();
  entity_info_ = entity_info;
  ProcessEnt(ent, compound_lib, mmcif_conform, entity_info_,
             atom_site_, pdbx_poly_seq_scheme_);
  this->Finalize(compound_lib);
}

void MMCifWriter::SetStructure(const ost::mol::EntityView& ent,
                               conop::CompoundLibPtr compound_lib,
                               bool mmcif_conform,
                               const std::vector<MMCifWriterEntity>& entity_info) {


  this->Setup();
  entity_info_ = entity_info;
  ProcessEnt(ent, compound_lib, mmcif_conform, entity_info_,
             atom_site_, pdbx_poly_seq_scheme_);
  this->Finalize(compound_lib);
}

void MMCifWriter::Setup() {
  if(structure_set_) {
    throw ost::io::IOException("SetStructure can be called only once on a "
                               "given MMCifWriter instance");
  }

  atom_type_ = Setup_atom_type_ptr();
  atom_site_ = Setup_atom_site_ptr();
  pdbx_poly_seq_scheme_ = Setup_pdbx_poly_seq_scheme_ptr();
  entity_ = Setup_entity_ptr();
  struct_asym_ = Setup_struct_asym_ptr();
  entity_poly_ = Setup_entity_poly_ptr();
  entity_poly_seq_ = Setup_entity_poly_seq_ptr();
  chem_comp_ = Setup_chem_comp_ptr();
  pdbx_entity_branch_ = Setup_pdbx_entity_branch_ptr();
}

void MMCifWriter::Finalize(ost::conop::CompoundLibPtr compound_lib) {
  // depending on the strategy (mmcif_conform), there might be gaps in the
  // entities mon_ids/ seq_olcs/ seq_can_olcs
  // The following function adds valid stuff depending on chain type
  // (e.g. UNK if we're having a peptide linking chain)
  ProcessUnknowns(entity_info_);

  Feed_entity(entity_, entity_info_);
  Feed_struct_asym(struct_asym_, entity_info_);
  Feed_entity_poly(entity_poly_, entity_info_);
  Feed_entity_poly_seq(entity_poly_seq_, entity_info_);
  Feed_chem_comp(chem_comp_, entity_info_, compound_lib);
  Feed_atom_type(atom_type_, atom_site_);
  Feed_pdbx_entity_branch(pdbx_entity_branch_, entity_info_);

  // finalize
  this->Push(entity_);
  this->Push(entity_poly_);
  this->Push(entity_poly_seq_);
  this->Push(chem_comp_);
  this->Push(struct_asym_);
  this->Push(pdbx_entity_branch_);
  this->Push(atom_type_);
  this->Push(atom_site_);
  this->Push(pdbx_poly_seq_scheme_);

  structure_set_ = true;
}

}} // ns
