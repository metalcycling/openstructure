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

#include <ost/mol/alg/biounit.hh>
#include <ost/mol/entity_view.hh>

namespace{

// dump and load vectors with various types of integers 
template<typename T>
void LoadIntVec(std::istream& stream, std::vector<T>& vec) {
  uint32_t size;
  stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  vec.resize(size);
  stream.read(reinterpret_cast<char*>(&vec[0]), size*sizeof(T));
}

template<typename T>
void DumpIntVec(std::ostream& stream, const std::vector<T>& vec) {
  uint32_t size = vec.size();
  stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  stream.write(reinterpret_cast<const char*>(&vec[0]), size*sizeof(T));
}

// dump and load strings
void Load(std::istream& stream, String& str) {
  uint32_t size;
  stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  str.resize(size);
  stream.read(&str[0], size);
}

void Dump(std::ostream& stream, const String& str) {
  uint32_t size = str.size();
  stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  stream.write(&str[0], size);
}

// dump and load vectors with strings
void Load(std::istream& stream, std::vector<String>& vec) {
  std::vector<uint8_t> string_sizes;
  String str;
  LoadIntVec(stream, string_sizes);
  Load(stream, str);
  vec.resize(string_sizes.size());
  int idx = 0;
  for(uint i = 0; i < string_sizes.size(); ++i) {
    vec[i] = str.substr(idx, string_sizes[i]);
    idx += string_sizes[i];
  }
}

void Dump(std::ostream& stream, const std::vector<String>& vec) {
  String total_string;
  std::vector<uint8_t> string_sizes;
  for(auto it = vec.begin(); it != vec.end(); ++it) {
    if(it->size() > std::numeric_limits<uint8_t>::max()) {
      std::stringstream ss;
      ss << "Max string size that can be encoded is "; 
      ss << std::numeric_limits<uint8_t>::max() << " cannot encode "<< *it; 
    }
    string_sizes.push_back(it->size());
    total_string += *it;
  }
  DumpIntVec(stream, string_sizes);
  Dump(stream, total_string);
}

// dump and load vectors with Mat4
void Load(std::istream& stream, std::vector<geom::Mat4>& vec) {
  uint32_t size;
  stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  vec.resize(size);
  for(uint i = 0; i < size; ++i) {
    stream.read(reinterpret_cast<char*>(vec[i].Data()),16*sizeof(Real));
  }
}

void Dump(std::ostream& stream, const std::vector<geom::Mat4>& vec) {
  uint32_t size = vec.size();
  stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  for(uint i = 0; i < size; ++i) {
    stream.write(reinterpret_cast<const char*>(vec[i].Data()), 16*sizeof(Real));
  }
}

} // anon ns


namespace ost{ namespace mol{ namespace alg{

BUInfo::BUInfo(const ost::io::MMCifInfoBioUnit& bu) {

  // translate MMCifInfoBioUnit objects into simpler data structures
  au_chains = bu.GetChainList();

  const std::vector<std::pair<int, int> >& bu_ch_intvl = bu.GetChainIntervalList();
  for(auto it = bu_ch_intvl.begin(); it != bu_ch_intvl.end(); ++it) {
    chain_intvl.push_back(it->first);      
    chain_intvl.push_back(it->second);      
  }

  const std::vector<std::vector<ost::io::MMCifInfoTransOpPtr> >& bu_op_list = bu.GetOperations();
  for(auto i = bu_op_list.begin(); i != bu_op_list.end(); ++i) {
    std::vector<geom::Mat4> mat_list;
    for(auto j = i->begin(); j != i->end(); ++j) {
      geom::Mat4 m;
      m.PasteRotation((*j)->GetMatrix());
      m.PasteTranslation((*j)->GetVector());
      mat_list.push_back(m);
    }
    operations.push_back(mat_list);
  }

  const std::vector<std::pair<int, int> >& bu_op_intvl = bu.GetOperationsIntervalList();
  for(auto it = bu_op_intvl.begin(); it != bu_op_intvl.end(); ++it) {
    op_intvl.push_back(it->first);      
    op_intvl.push_back(it->second);      
  }
}

void BUInfo::ToStream(std::ostream& stream) const {
  Dump(stream, au_chains);
  DumpIntVec(stream, chain_intvl);
  uint32_t size = operations.size();
  stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  for(auto it = operations.begin(); it != operations.end(); ++it) {
    Dump(stream, *it);
  }
  DumpIntVec(stream, op_intvl);
}

BUInfo BUInfo::FromStream(std::istream& stream) {
  BUInfo info;
  Load(stream, info.au_chains);
  LoadIntVec(stream, info.chain_intvl);
  uint32_t size = 0;
  stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  info.operations.resize(size);
  for(uint i = 0; i < size; ++i) {
    Load(stream, info.operations[i]);
  }
  LoadIntVec(stream, info.op_intvl);
  return info;
}

BUInfo BUInfo::FromString(const String& s) {
  std::istringstream in_stream(s);
  BUInfo info = BUInfo::FromStream(in_stream);
  return info;
}

String BUInfo::ToString() const {
  std::ostringstream out_stream;
  this->ToStream(out_stream);
  return out_stream.str();
}

const std::vector<std::vector<String> >& BUInfo::GetAUChains() const {
  if(au_chains_.empty()) {
    this->_InitTransforms();
  }
  return au_chains_;
}

const std::vector<std::vector<geom::Mat4> >& BUInfo::GetTransformations() const {
  if(transforms_.empty()) {
    this->_InitTransforms();
  }
  return transforms_;
}

void BUInfo::_InitTransforms() const {
  int n_intervals = chain_intvl.size() / 2;
  for(int intvl_idx = 0; intvl_idx < n_intervals; ++intvl_idx) {
    // extract relevant chain names from asu
    ////////////////////////////////////////
    std::vector<String> chain_names;
    int chain_start = chain_intvl[2*intvl_idx];
    int chain_end = chain_intvl[2*intvl_idx+1];
    for(int ch_idx = chain_start; ch_idx < chain_end; ++ch_idx) {
      chain_names.push_back(au_chains[ch_idx]);
    }
    au_chains_.push_back(chain_names);
    // extract operations that will be applied to those chains
    //////////////////////////////////////////////////////////
    std::vector<geom::Mat4> transforms;
    int op_start = op_intvl[2*intvl_idx];
    int op_end = op_intvl[2*intvl_idx+1];
    if(op_end > op_start) {
      for(auto it = operations[op_start].begin();
          it != operations[op_start].end(); ++it) {
        transforms.push_back(*it);
      }
      ++op_start;
      while(op_start < op_end) {
        std::vector<geom::Mat4> tmp_transforms;
        for(auto i = operations[op_start].begin(); 
            i != operations[op_start].end(); ++i) {
          for(auto j = transforms.begin(); j != transforms.end(); ++j) {
            tmp_transforms.push_back((*j)*(*i));
          }
        }
        transforms = tmp_transforms;
        ++op_start;
      }
    }
    transforms_.push_back(transforms);
  }
}

ost::mol::EntityHandle CreateBU(const ost::mol::EntityHandle& asu,
                                const ost::io::MMCifInfoBioUnit& bu) {
  BUInfo bu_info(bu);
  return CreateBU(asu, bu_info);
}

ost::mol::EntityHandle CreateBU(const ost::mol::EntityHandle& asu,
                                const BUInfo& bu_info) {
  ost::mol::EntityHandle ent = ost::mol::CreateEntity();
  ent.SetName(asu.GetName());
  ost::mol::XCSEditor ed = ent.EditXCS(mol::BUFFERED_EDIT);

  // For chain naming. First copy with transformation: 2.<au_cname>, second
  // 3.<au_cname> etc.
  std::map<String, int> chain_counter;

  // The name 1.<au_cname> is reserved for that particular AU chain with
  // identity transform, i.e. the copy of the actual AU chain. We need to keep
  // track of this as there can only be one.
  std::set<String> au_chain_copies;

  const std::vector<std::vector<String> >& au_chains = bu_info.GetAUChains();
  const std::vector<std::vector<geom::Mat4> >& transforms =
  bu_info.GetTransformations();

  for(uint chain_intvl = 0; chain_intvl < au_chains.size(); ++chain_intvl) {
    if(au_chains[chain_intvl].empty()) continue;
    // derive all bonds related to that chain_intvl
    // potentially also interchain bonds
    std::stringstream query_ss;
    query_ss << "cname=" << au_chains[chain_intvl][0];
    for(uint i = 1; i < au_chains[chain_intvl].size(); ++i) {
      query_ss << ',' << au_chains[chain_intvl][i];
    }
    ost::mol::EntityView asu_view = asu.Select(query_ss.str());
    const ost::mol::BondHandleList& bond_list = asu_view.GetBondList();

    // process all transformations
    for(uint t_idx = 0; t_idx < transforms[chain_intvl].size(); ++t_idx) {
      const geom::Mat4& m = transforms[chain_intvl][t_idx];
      // check if m is identity matrix => no transformation applied
      bool is_identity = true;
      geom::Mat4 identity_matrix = geom::Mat4::Identity();
      const Real* m_data = m.Data();
      const Real* identity_data = identity_matrix.Data();
      for(int i = 0; i < 16; ++i) {
        if(std::abs(m_data[i] - identity_data[i]) > 1e-5) {
          is_identity = false;
          break;
        }
      }

      // key: au_at.GetHashCode, value: bu_at
      // required for bond buildup in the end
      std::map<long, AtomHandle> atom_mapper;
      for(uint c_idx = 0; c_idx < au_chains[chain_intvl].size(); ++c_idx) {
        String au_cname = au_chains[chain_intvl][c_idx];

        std::stringstream bu_cname_ss;
        if(is_identity && au_chain_copies.find(au_cname) == au_chain_copies.end()) {
          bu_cname_ss << "1." << au_cname; // 1.<au_cname> reserved for AU chain
                                           // without transformation
                                           // at least the first of it...
                                           // as of January 2024, there were 3
                                           // entries (8qn6, 8x1h, 2c0x) where
                                           // the identity transform is applied
                                           // more than once on the same AU
                                           // chain, effectively leading to
                                           // chains sitting on top of each
                                           // other... But hey, bullshit in,
                                           // bullshit out
          au_chain_copies.insert(au_cname);
        } else {
          if(chain_counter.find(au_cname) == chain_counter.end()) {
            chain_counter[au_cname] = 2;
          }
          bu_cname_ss << chain_counter[au_cname] << '.' << au_cname;
          chain_counter[au_cname] += 1;
        }
        ost::mol::ChainHandle asu_ch = asu.FindChain(au_cname);
        if(!asu_ch.IsValid()) {
          std::stringstream ss;
          ss << "Cannot construct biounit with asu chain "<<au_cname;
          ss << ". Specified interval only has: " <<au_chains[chain_intvl][0];
          for(uint i = 1; i < au_chains[chain_intvl].size(); ++i) {
            ss << ',' << au_chains[chain_intvl][i];
          }
          throw ost::Error(ss.str());
        }
        ost::mol::ChainHandle bu_ch = ed.InsertChain(bu_cname_ss.str());
        ed.SetChainType(bu_ch, asu_ch.GetType());
        ost::mol::ResidueHandleList au_res_list = asu_ch.GetResidueList();
        for(auto res_it = au_res_list.begin();
            res_it != au_res_list.end(); ++res_it) {
          ost::mol::ResidueHandle bu_res = ed.AppendResidue(bu_ch,
            res_it->GetName(), res_it->GetNumber());
          bu_res.SetOneLetterCode(res_it->GetOneLetterCode());
          bu_res.SetSecStructure(res_it->GetSecStructure());
          bu_res.SetChemClass(res_it->GetChemClass());
          bu_res.SetChemType(res_it->GetChemType());
          bu_res.SetIsProtein(res_it->IsProtein());
          bu_res.SetIsLigand(res_it->IsLigand());
          ost::mol::AtomHandleList au_at_list = res_it->GetAtomList();
          for(auto at_it = au_at_list.begin(); at_it != au_at_list.end(); ++at_it) {
            geom::Vec3 bu_at_pos = geom::Vec3(m*geom::Vec4(at_it->GetPos()));
            ost::mol::AtomHandle bu_at = ed.InsertAtom(bu_res, at_it->GetName(),
                                                       bu_at_pos,
                                                       at_it->GetElement(),
                                                       at_it->GetOccupancy(),
                                                       at_it->GetBFactor(),
                                                       at_it->IsHetAtom());
            atom_mapper[at_it->GetHashCode()] = bu_at;
          }
        }
      }

      // connect
      for(auto it = bond_list.begin(); it != bond_list.end(); ++it) {
        ed.Connect(atom_mapper[it->GetFirst().GetHashCode()],
                   atom_mapper[it->GetSecond().GetHashCode()]);
      }

    }
  }
  return ent;
}

}}} // ns
