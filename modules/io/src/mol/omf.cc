#include <algorithm>

#include <ost/mol/atom_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/xcs_editor.hh>

#include "omf.hh"

namespace{

  void ConstructOPos(const geom::Vec3& ca_pos, const geom::Vec3& c_pos,
                     const geom::Vec3& n_next_pos, geom::Vec3& o_pos) {
    geom::Vec3 o_vec = geom::Normalize(geom::Normalize(c_pos - ca_pos) +
                                       geom::Normalize(c_pos - n_next_pos));
    o_pos = c_pos + 1.230*o_vec;
  }

  void ConstructAtomPos(const geom::Vec3& A, const geom::Vec3& B,
                        const geom::Vec3& C, Real bond_length, Real angle,
                        Real dihedral, geom::Vec3& pos) {

    Real x,xx,x1;
    x = std::tan((Real(M_PI)-angle)*Real(0.5));
    xx = x*x;
    x1 = Real(1.0)+xx;
    Real bond_length_x_sin_angle = bond_length*Real(2.0)*x/x1;
    Real bond_length_x_cos_angle = bond_length*(Real(1.0)-xx)/x1;

    x = std::tan(Real(0.5)*dihedral);
    xx = x*x;
    x1 = Real(1.0)+xx;
    Real sin_dihedral = Real(2.0)*x/x1;
    Real cos_dihedral = (Real(1.0)-xx)/x1;

    Real ab[] = {B[0]-A[0],B[1]-A[1],B[2]-A[2]};

    Real norm_bc[] = {C[0]-B[0],C[1]-B[1],C[2]-B[2]};
    Real a = norm_bc[0] * norm_bc[0];
    a += norm_bc[1] * norm_bc[1];
    a += norm_bc[2] * norm_bc[2];
    a = Real(1.0) / std::sqrt(a);

    norm_bc[0]*=a;
    norm_bc[1]*=a;
    norm_bc[2]*=a;

    Real norm_n[] = {ab[1]*norm_bc[2]-norm_bc[1]*ab[2],
                     ab[2]*norm_bc[0]-norm_bc[2]*ab[0],
                     ab[0]*norm_bc[1]-norm_bc[0]*ab[1]};

    a = norm_n[0]*norm_n[0];  
    a += norm_n[1]*norm_n[1];  
    a += norm_n[2]*norm_n[2];  
    a = Real(1.0) / std::sqrt(a);

    norm_n[0]*=a;
    norm_n[1]*=a;
    norm_n[2]*=a;

    Real n_x_bc[] = {norm_n[1]*norm_bc[2]-norm_bc[1]*norm_n[2],
                     norm_n[2]*norm_bc[0]-norm_bc[2]*norm_n[0],
                     norm_n[0]*norm_bc[1]-norm_bc[0]*norm_n[1]};

    pos[0] = norm_bc[0]*bond_length_x_cos_angle + 
             n_x_bc[0]*cos_dihedral*bond_length_x_sin_angle + 
             norm_n[0]*sin_dihedral*bond_length_x_sin_angle + C[0];

    pos[1] = norm_bc[1]*bond_length_x_cos_angle + 
             n_x_bc[1]*cos_dihedral*bond_length_x_sin_angle + 
             norm_n[1]*sin_dihedral*bond_length_x_sin_angle + C[1];

    pos[2] = norm_bc[2]*bond_length_x_cos_angle + 
             n_x_bc[2]*cos_dihedral*bond_length_x_sin_angle + 
             norm_n[2]*sin_dihedral*bond_length_x_sin_angle + C[2];
  }

  // some hash function we need for an unordered_map
  // stolen from https://stackoverflow.com/questions/32685540/why-cant-i-compile-an-unordered-map-with-a-pair-as-key
  struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
      auto h1 = std::hash<T1>{}(p.first);
      auto h2 = std::hash<T2>{}(p.second);
      // Mainly for demonstration purposes, i.e. works but is overly simple
      // In the real world, use sth. like boost.hash_combine
      return h1 ^ h2;  
    }
  };

  // define hash function, so we can use ResidueDefinition as key in an
  // unordered map. The used hash function is overly simple and gives a hash 
  // collision whenever we have two residues of same name but different atom
  // composition. That's hopefully rare...
  struct ResidueDefinitionHash {
    std::size_t operator()(const ost::io::ResidueDefinition& r) const {
      return std::hash<String>()(r.name);
    }
  };

  // some helpers
  void RealToIntVec(const std::vector<Real>& real_vec, 
                    std::vector<int>& int_vec, Real factor) {
    int_vec.resize(real_vec.size());
    for(uint i = 0; i < real_vec.size(); ++i) {
      int_vec[i] = std::round(factor*real_vec[i]);
    }
  }

  void IntToRealVec(const std::vector<int>& int_vec, 
                    std::vector<Real>& real_vec, Real factor) {
    real_vec.resize(int_vec.size());
    for(uint i = 0; i < int_vec.size(); ++i) {
      real_vec[i] = factor*int_vec[i];
    }
  }

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

  // delta/runlength encodings/decodings
  void DeltaEncoding(const std::vector<int>& in, std::vector<int>& out) {
    out.clear();
    if(in.empty()) {
      return;
    }
    out.reserve(in.size());
    out.push_back(in[0]);
    for(uint i = 1; i < in.size(); ++i) {
      out.push_back(in[i] - in[i-1]);
    }
  }

  void DeltaDecoding(const std::vector<int>& in, std::vector<int>& out) {
    out.clear();
    if(in.empty()) {
      return;
    }
    out.reserve(in.size());
    out.push_back(in[0]);
    for(uint i = 1; i < in.size(); ++i) {
      out.push_back(out.back() + in[i]);
    }
  }

  void RunLengthEncoding(const std::vector<int>& in, std::vector<int>& out) {
    out.clear();
    if(in.empty()) {
      return;
    }
    int current_item = in[0];
    int run_length = 1;
    for(uint i = 1; i < in.size(); ++i) {
      if(in[i] == current_item) {
        ++run_length;
      } else {
        out.push_back(current_item);
        out.push_back(run_length);
        current_item = in[i];
        run_length = 1;
      }
    }
    out.push_back(current_item);
    out.push_back(run_length);
  }

  void RunLengthDecoding(const std::vector<int>& in, std::vector<int>& out) {
    out.clear();
    if(in.empty()) {
      return;
    }
    int n_runs = in.size() / 2;
    for(int run_idx = 0; run_idx < n_runs; ++run_idx) {
      int value = in[run_idx*2];
      int run_length = in[run_idx*2+1];
      out.insert(out.end(), run_length, value);
    }  
  }

  // functionality to perform integer packing
  int IntegerPackingSize(const std::vector<int>& vec, int min, int max) {
    // This function only estimates the number of elements in integer
    // packing without actually doing it.

    if(max <= min) {
      throw ost::Error("Min max error in IntegerPackingSize");
    }

    // We don't allow unsigned packing here => min must be negative,
    // max must be positive
    if(min >= 0) {
      throw ost::Error("Min val in IntegerPacking must be negative");
    }

    if(max <= 0) {
      throw ost::Error("Max val in IntegerPacking must be positive");
    }

    int n = 0;
    int abs_min = std::abs(min);
    for(auto it = vec.begin(); it != vec.end(); ++it) {
      if(*it >= max) {
        n += *it/max;
        ++n; 
      } else if (*it <= min) {
        n += std::abs(*it)/abs_min;
        ++n;
      } else {
        ++n;
      }
    }
    return n;
  }

  template<typename T>
  void IntegerPacking(const std::vector<int>& in, std::vector<T>& out,
                      int min, int max) {

    if(min < std::numeric_limits<T>::min()) {
      throw ost::Error("Invalid min val in IntegerPacking");
    }

    if(max > std::numeric_limits<T>::max()) {
      throw ost::Error("Invalid max val in IntegerPacking");
    }

    if(max <= min) {
      throw ost::Error("Min max error in IntegerPackingSize");
    }

    // We don't allow unsigned packing here => min must be negative,
    // max must be positive
    if(min >= 0) {
      throw ost::Error("Min val in IntegerPacking must be negative");
    }

    if(max <= 0) {
      throw ost::Error("Max val in IntegerPacking must be positive");
    }

    out.clear();
    int abs_min = std::abs(min);
    for(auto it = in.begin(); it != in.end(); ++it) {
      if(*it >= max) {
        int n = *it/max;
        out.insert(out.end(), n, max);
        out.push_back(*it - n*max); 
      } else if (*it <= min) {
        int n = std::abs(*it)/abs_min;
        out.insert(out.end(), n, min);
        out.push_back(*it + n*abs_min);
      } else {
        out.push_back(*it);
      }
    }
  }

  template<typename T>
  void IntegerUnpacking(const std::vector<T>& in, std::vector<int>& out,
                        int min, int max) {

    if(min < std::numeric_limits<T>::min()) {
      throw ost::Error("Invalid min val in IntegerUnpacking");
    }

    if(max > std::numeric_limits<T>::max()) {
      throw ost::Error("Invalid max val in IntegerUnpacking");
    }

    if(max <= min) {
      throw ost::Error("Min max error in IntegerUnpacking");
    }

    // We don't allow unsigned packing here => min must be negative,
    // max must be positive
    if(min >= 0) {
      throw ost::Error("Min val in IntegerUnpacking must be negative");
    }

    if(max <= 0) {
      throw ost::Error("Max val in IntegerUnpacking must be positive");
    }

    out.clear();
    int idx = 0;
    int in_size = in.size();
    while(idx < in_size) {
      int val = in[idx];
      if(val == max) {
        int summed_val = val;
        ++idx;
        while(true) {
          summed_val += in[idx];
          if(in[idx] != max) {
            ++idx;
            break;
          }
          ++idx;
        }
        out.push_back(summed_val);
      } else if (val == min) {
        int summed_val = min;
        ++idx;
        while(true) {
          summed_val += in[idx];
          if(in[idx] != min) {
            ++idx;
            break;
          }
          ++idx;
        }
        out.push_back(summed_val);
      } else {
        out.push_back(val);
        ++idx;
      }
    }
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

  // dump and load chars
  void Load(std::istream& stream, char& ch) {
    stream.read(&ch, sizeof(char));
  }

  void Dump(std::ostream& stream, char ch) {
    stream.write(&ch, sizeof(char));
  }

  // dump and load maps with string as key and ChainDataPtr as value
  void Load(std::istream& stream, 
            std::map<String, ost::io::ChainDataPtr>& map,
            const std::vector<ost::io::ResidueDefinition>& res_def,
            int version, bool lossy, bool avg_bfactors, bool round_bfactors,
            bool skip_ss, bool infer_aa_pos) {
    uint32_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    map.clear();
    for(uint i = 0; i < size; ++i) {
      ost::io::ChainDataPtr p(new ost::io::ChainData);
      p->FromStream(stream, res_def, version, lossy, avg_bfactors,
                    round_bfactors, skip_ss, infer_aa_pos);
      map[p->ch_name] = p;
    }
  }

  void Dump(std::ostream& stream, 
            const std::map<String, ost::io::ChainDataPtr>& map,
            const std::vector<ost::io::ResidueDefinition>& res_def,
            bool lossy, bool avg_bfactors, bool round_bfactors,
            bool skip_ss, bool infer_aa_pos) {
    uint32_t size = map.size();
    stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    for(auto it = map.begin(); it != map.end(); ++it) {
        // we don't dump the key (chain name), that's an attribute of the
        // chain itself anyway
      it->second->ToStream(stream, res_def, lossy, avg_bfactors,
                           round_bfactors, skip_ss, infer_aa_pos); 
    }
  }

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

  void Load(std::istream& stream, std::vector<int>& vec) {
    int8_t encoding;
    stream.read(reinterpret_cast<char*>(&encoding), sizeof(int8_t));
    if(encoding == 8) {
      std::vector<int8_t> int8_vec;
      LoadIntVec(stream, int8_vec);
      IntegerUnpacking(int8_vec, vec, std::numeric_limits<int8_t>::min(),
                       std::numeric_limits<int8_t>::max());
    } else if(encoding == 16) {
      std::vector<int16_t> int16_vec;
      LoadIntVec(stream, int16_vec);
      IntegerUnpacking(int16_vec, vec, std::numeric_limits<int16_t>::min(),
                       std::numeric_limits<int16_t>::max());
    } else if(encoding == 32) {
      LoadIntVec(stream, vec);
    } else {
      throw ost::Error("Encountered unknown encoding when loading int vec." );
    }
  }

  void Dump(std::ostream& stream, const std::vector<int>& vec) {

    int8_t encoding = 32;

    // check whether we can pack tighter
    int n_16 = IntegerPackingSize(vec, std::numeric_limits<int16_t>::min(), 
                                  std::numeric_limits<int16_t>::max());
    if(n_16*sizeof(int16_t) < vec.size()*sizeof(int)) {
      // less bytes required...
      encoding = 16;
      //even tighter?
      int n_8 = IntegerPackingSize(vec, std::numeric_limits<int8_t>::min(), 
                                   std::numeric_limits<int8_t>::max());
      if(n_8*sizeof(int8_t) < n_16*sizeof(int16_t)) {
        encoding = 8;
      }
    }

    stream.write(reinterpret_cast<char*>(&encoding), sizeof(int8_t));

    if(encoding == 32) {
      DumpIntVec(stream, vec);
    } else if(encoding == 16) {
      std::vector<int16_t> packed;
      IntegerPacking(vec, packed, std::numeric_limits<int16_t>::min(),
                     std::numeric_limits<int16_t>::max());
      DumpIntVec(stream, packed);
    } else if(encoding == 8) {
      std::vector<int8_t> packed;
      IntegerPacking(vec, packed, std::numeric_limits<int8_t>::min(),
                     std::numeric_limits<int8_t>::max());
      DumpIntVec(stream, packed);
    } else {
      throw ost::Error("AAAAAAAAaaaaaa!");
    }
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

  // dump and load vectors with ResidueDefinition
  void Load(std::istream& stream, std::vector<ost::io::ResidueDefinition>& vec) {
    uint32_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    vec.resize(size);
    for(uint i = 0; i < size; ++i) {
      vec[i].FromStream(stream);
    }
  }

  void Dump(std::ostream& stream,
            const std::vector<ost::io::ResidueDefinition>& vec) {
    uint32_t size = vec.size();
    stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    for(uint i = 0; i < size; ++i) {
      vec[i].ToStream(stream);
    }
  }

  // dump and load vectors with BioUnitDefinition
  void Load(std::istream& stream,
            std::vector<ost::io::BioUnitDefinition>& vec) {
    uint32_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    vec.resize(size);
    for(uint i = 0; i < size; ++i) {
      vec[i].FromStream(stream);
    }
  }

  void Dump(std::ostream& stream,
            const std::vector<ost::io::BioUnitDefinition>& vec) {
    uint32_t size = vec.size();
    stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    for(uint i = 0; i < size; ++i) {
      vec[i].ToStream(stream);
    }
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
      stream.write(reinterpret_cast<const char*>(vec[i].Data()),
                   16*sizeof(Real));
    }
  }

  // knowledge based dump and load functions, i.e. awesome compressions
  void LoadIsHetatm(std::istream& stream, std::vector<bool>& vec) {
    std::vector<int> run_length_encoded;
    Load(stream, run_length_encoded);
    std::vector<int> int_vec;
    RunLengthDecoding(run_length_encoded, int_vec);
    vec.assign(int_vec.begin(), int_vec.end());
  }

  void DumpIsHetatm(std::ostream& stream, const std::vector<bool>& vec) {
    std::vector<int> int_vec(vec.begin(), vec.end());
    std::vector<int> run_length_encoded;
    RunLengthEncoding(int_vec, run_length_encoded);
    Dump(stream, run_length_encoded);
  }

  void LoadBonds(std::istream& stream, std::vector<int>& vec) {
    std::vector<int> delta_encoded;
    Load(stream, delta_encoded);
    DeltaDecoding(delta_encoded, vec);
  }

  void DumpBonds(std::ostream& stream, const std::vector<int>& vec) {
    std::vector<int> delta_encoded;
    DeltaEncoding(vec, delta_encoded);
    Dump(stream, delta_encoded);
  }

  void LoadBondOrders(std::istream& stream, std::vector<int>& vec) {
    std::vector<int> run_length_encoded;
    Load(stream, run_length_encoded);
    RunLengthDecoding(run_length_encoded, vec);
  }

  void DumpBondOrders(std::ostream& stream, const std::vector<int>& vec) {
    std::vector<int> run_length_encoded;
    RunLengthEncoding(vec, run_length_encoded);
    Dump(stream, run_length_encoded);
  }

  void LoadPosVec(std::istream& stream, std::vector<Real>& vec, bool lossy) {
    std::vector<int> delta_encoded;
    Load(stream, delta_encoded);
    std::vector<int> int_vec;
    DeltaDecoding(delta_encoded, int_vec);
    if(lossy) {
      IntToRealVec(int_vec, vec, 0.1);  
    } else {
      IntToRealVec(int_vec, vec, 0.001);  
    }
  }

  void LoadPositions(std::istream& stream, geom::Vec3List& positions,
                     bool lossy) {
    std::vector<Real> x_pos;
    std::vector<Real> y_pos;
    std::vector<Real> z_pos;
    LoadPosVec(stream, x_pos, lossy);
    LoadPosVec(stream, y_pos, lossy);
    LoadPosVec(stream, z_pos, lossy);
    positions.resize(x_pos.size());
    for(uint i = 0; i < positions.size(); ++i) {
      positions[i] = geom::Vec3(x_pos[i], y_pos[i], z_pos[i]);
    }
  }

  void DumpPosVec(std::ostream& stream, const std::vector<Real>& vec,
                  bool lossy) {
    std::vector<int> int_vec;
    if(lossy) {
      RealToIntVec(vec, int_vec, 10);  
    } else {
      RealToIntVec(vec, int_vec, 1000);
    }
    std::vector<int> delta_compressed;
    DeltaEncoding(int_vec, delta_compressed);
    Dump(stream, delta_compressed);    
  }

  void DumpPositions(std::ostream& stream, const geom::Vec3List& positions,
                     bool lossy) {
    std::vector<Real> x_pos(positions.size());
    std::vector<Real> y_pos(positions.size());
    std::vector<Real> z_pos(positions.size());
    for(uint i = 0; i < positions.size(); ++i) {
      x_pos[i] = positions[i][0];
      y_pos[i] = positions[i][1];
      z_pos[i] = positions[i][2];
    }
    DumpPosVec(stream, x_pos, lossy);
    DumpPosVec(stream, y_pos, lossy);
    DumpPosVec(stream, z_pos, lossy);
  }

  void DumpDihedrals(std::ostream& stream, const std::vector<Real>& dihedrals) {
    std::vector<int> int_vec(dihedrals.size());
    for(size_t i = 0; i < dihedrals.size(); ++i) {
      int_vec[i] = round((dihedrals[i] + M_PI)/(2*M_PI)*255);
    }
    Dump(stream, int_vec);
  }

  void LoadDihedrals(std::istream& stream, std::vector<Real>& dihedrals) {
    std::vector<int> int_vec;
    Load(stream, int_vec);
    dihedrals.resize(int_vec.size());
    for(size_t i = 0; i < int_vec.size(); ++i) {
      dihedrals[i] = static_cast<Real>(int_vec[i])/255*2*M_PI-M_PI;
    }
  }

  void LoadBFactors(std::istream& stream, std::vector<Real>& bfactors,
                    bool round_bfactors) {

    int8_t bfactor_encoding = 0;
    stream.read(reinterpret_cast<char*>(&bfactor_encoding), sizeof(int8_t));
    if(bfactor_encoding == 0) {
      std::vector<int> delta_encoded;
      Load(stream, delta_encoded);
      std::vector<int> int_vec;
      DeltaDecoding(delta_encoded, int_vec);
      if(round_bfactors) {
        IntToRealVec(int_vec, bfactors, 1.0);
      } else {
        IntToRealVec(int_vec, bfactors, 0.01);
      }
    } else if(bfactor_encoding == 42) {
      std::vector<int> runlength_encoded;
      Load(stream, runlength_encoded);
      std::vector<int> int_vec;
      RunLengthDecoding(runlength_encoded, int_vec);
      if(round_bfactors) {
        IntToRealVec(int_vec, bfactors, 1.0);
      } else {
        IntToRealVec(int_vec, bfactors, 0.01);
      }
    } else {
      throw ost::Error("Observed invalid bfactor encoding");
    }
  }

  void DumpBFactors(std::ostream& stream, const std::vector<Real>& bfactors,
                    bool round_bfactors) {
    std::vector<int> int_vec;
    if(round_bfactors) {
      RealToIntVec(bfactors, int_vec, 1);
    } else {
      RealToIntVec(bfactors, int_vec, 100);
    }

    // Hack: some structures (e.g. EM) have all bfactors set to 0.0
    // this efficiently compresses with runlength encoding.
    // Let's sacrifice a byte to mark that
    std::vector<int> run_length_encoded;
    RunLengthEncoding(int_vec, run_length_encoded);
    if(static_cast<float>(run_length_encoded.size())/int_vec.size() < 0.42) {
      int8_t bfactor_encoding = 42;
      stream.write(reinterpret_cast<char*>(&bfactor_encoding), sizeof(int8_t));
      Dump(stream, run_length_encoded);
    } else {
      // continue with delta encoding
      int8_t bfactor_encoding = 0;
      stream.write(reinterpret_cast<char*>(&bfactor_encoding), sizeof(int8_t));
      std::vector<int> delta_encoded;
      DeltaEncoding(int_vec, delta_encoded);
      Dump(stream, delta_encoded);
    }
  }

  void LoadOccupancies(std::istream& stream, std::vector<Real>& occupancies) {
    std::vector<int> run_length_encoded;
    Load(stream, run_length_encoded);
    std::vector<int> int_vec;
    RunLengthDecoding(run_length_encoded, int_vec);
    IntToRealVec(int_vec, occupancies, 0.01);
  }

  void DumpOccupancies(std::ostream& stream,
                       const std::vector<Real>& occupancies) {
    std::vector<int> int_vec;
    RealToIntVec(occupancies, int_vec, 100.0);
    std::vector<int> run_length_encoded;
    RunLengthEncoding(int_vec, run_length_encoded);
    Dump(stream, run_length_encoded);
  }

  void LoadInsertionCodes(std::istream& stream, std::vector<char>& ins_codes) {
    std::vector<int> run_length_encoded;
    Load(stream, run_length_encoded);
    std::vector<int> int_vec;
    RunLengthDecoding(run_length_encoded, int_vec);
    ins_codes.assign(int_vec.begin(), int_vec.end());
  }

  void DumpInsertionCodes(std::ostream& stream,
                          const std::vector<char>& ins_codes) {
    std::vector<int> int_vec(ins_codes.begin(), ins_codes.end());
    std::vector<int> run_length_encoded;
    RunLengthEncoding(int_vec, run_length_encoded);
    Dump(stream, run_length_encoded);
  }

  void LoadRnums(std::istream& stream, std::vector<int>& rnums) {
    std::vector<int> run_length_encoded;
    Load(stream, run_length_encoded);
    std::vector<int> delta_encoded;
    RunLengthDecoding(run_length_encoded, delta_encoded);
    DeltaDecoding(delta_encoded, rnums);
  }

  void DumpRnums(std::ostream& stream, const std::vector<int>& rnums) {
    std::vector<int> delta_encoded;
    DeltaEncoding(rnums, delta_encoded);
    std::vector<int> run_length_encoded;
    RunLengthEncoding(delta_encoded, run_length_encoded);
    Dump(stream, run_length_encoded);
  }

  void LoadResDefIndices(std::istream& stream, std::vector<int>& indices) {
    Load(stream, indices);
  }

  void DumpResDefIndices(std::ostream& stream,
                         const std::vector<int>& indices) {
    Dump(stream, indices);
  }

  void LoadSecStructures(std::istream& stream,
                         std::vector<char>& sec_structures) {
    std::vector<int> run_length_encoded;
    Load(stream, run_length_encoded);
    std::vector<int> transformed_sec_structures;
    RunLengthDecoding(run_length_encoded, transformed_sec_structures);
    sec_structures.clear();
    sec_structures.reserve(transformed_sec_structures.size());
    for(auto it = transformed_sec_structures.begin();
        it != transformed_sec_structures.end(); ++it) {
      switch(*it) {
        case 0: {
          sec_structures.push_back(ost::mol::SecStructure::ALPHA_HELIX);
          break;
        }
        case 1: {
          sec_structures.push_back(ost::mol::SecStructure::COIL);
          break;
        }
        case 2: {
          sec_structures.push_back(ost::mol::SecStructure::THREE_TEN_HELIX);
          break;
        }
        case 3: {
          sec_structures.push_back(ost::mol::SecStructure::TURN);
          break;
        }
        case 4: {
          sec_structures.push_back(ost::mol::SecStructure::EXTENDED);
          break;
        }
        case 5: {
          sec_structures.push_back(ost::mol::SecStructure::BETA_BRIDGE);
          break;
        }
        case 6: {
          sec_structures.push_back(ost::mol::SecStructure::BEND);
          break;
        }
        case 7: {
          sec_structures.push_back(ost::mol::SecStructure::PI_HELIX);
          break;
        }
        default: {
          throw ost::Error("Invalid sec structure observed");
        }
      }
    }
  }

  void DumpSecStructures(std::ostream& stream,
                         const std::vector<char>& sec_structures) {
    std::vector<int> transformed_sec_structures;
    transformed_sec_structures.reserve(sec_structures.size());
    for(auto it = sec_structures.begin(); it != sec_structures.end(); ++it) {
      switch(*it) {
        case ost::mol::SecStructure::ALPHA_HELIX: {
          transformed_sec_structures.push_back(0);
          break;
        }
        case ost::mol::SecStructure::COIL: {
          transformed_sec_structures.push_back(1);
          break;
        }
        case ost::mol::SecStructure::THREE_TEN_HELIX: {
          transformed_sec_structures.push_back(2);
          break;
        }
        case ost::mol::SecStructure::TURN: {
          transformed_sec_structures.push_back(3);
          break;
        }
        case ost::mol::SecStructure::EXTENDED: {
          transformed_sec_structures.push_back(4);
          break;
        }
        case ost::mol::SecStructure::BETA_BRIDGE: {
          transformed_sec_structures.push_back(5);
          break;
        }
        case ost::mol::SecStructure::BEND: {
          transformed_sec_structures.push_back(6);
          break;
        }
        case ost::mol::SecStructure::PI_HELIX: {
          transformed_sec_structures.push_back(7);
          break;
        }
        default: {
          throw ost::Error("Invalid sec structure observed");
        }
      }
    }
    std::vector<int> run_length_encoded;
    RunLengthEncoding(transformed_sec_structures, run_length_encoded);
    Dump(stream, run_length_encoded);
  }

  void DumpName(std::ostream& stream, const String& name) {
    if(name.size() > std::numeric_limits<uint8_t>::max()) {
      std::stringstream ss;
      ss << "Max name size that can be dumped is ";
      ss << std::numeric_limits<uint8_t>::max << ". ";
      ss << "got: "<<name<<std::endl;
      throw ost::Error(ss.str());
    }
    uint8_t size = name.size();
    stream.write(reinterpret_cast<char*>(&size), sizeof(uint8_t));
    stream.write(reinterpret_cast<const char*>(&name[0]), size*sizeof(char));
  }

  void LoadName(std::istream& stream, String& name) {
    uint8_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(uint8_t));
    name.resize(size);
    stream.read(reinterpret_cast<char*>(&name[0]), size*sizeof(uint8_t));
  }
}


namespace ost { namespace io {

ResidueDefinition::ResidueDefinition(const ost::mol::ResidueHandle& res) {
  name = res.GetName();
  olc = res.GetOneLetterCode();
  chem_type = char(res.GetChemType());
  chem_class = char(res.GetChemClass());
  rotamer_setup = false;
  ost::mol::AtomHandleList at_list = res.GetAtomList();
  for(auto it = at_list.begin(); it != at_list.end(); ++it) {
    anames.push_back(it->GetName());
  }

  // sort atom names and store other info accordingly
  // sorting is required to compare definitions from different residues
  std::sort(anames.begin(), anames.end());
  std::unordered_map<unsigned long, int> hash_idx_mapper;
  for(auto it = anames.begin(); it != anames.end(); ++it) {
    ost::mol::AtomHandle at = res.FindAtom(*it);
    hash_idx_mapper[at.GetHashCode()] = elements.size();
    elements.push_back(at.GetElement());
    is_hetatm.push_back(at.IsHetAtom());
  }

  // in principle we want to store all unique bonds, but nevertheless
  // use a map here to track the bond orders at the same time
  std::unordered_map<std::pair<int, int>,  int, pair_hash> bond_order_map;
  for(auto at_it = at_list.begin(); at_it != at_list.end(); ++at_it) {
    ost::mol::BondHandleList bond_list = at_it->GetBondList();  
    for(auto bond_it = bond_list.begin(); bond_it != bond_list.end();
        ++bond_it) {
      // The atom represented by at_it is either first OR second atom of that
      // bond. So we check whether BOTH are in hash_idx_mapper to exclude bonds
      // to other residues.
      unsigned long h1 = bond_it->GetFirst().GetHashCode();
      unsigned long h2 = bond_it->GetSecond().GetHashCode();
      if(hash_idx_mapper.find(h1) != hash_idx_mapper.end() &&
         hash_idx_mapper.find(h2) != hash_idx_mapper.end()) {
        int i1 = hash_idx_mapper[h1];
        int i2 = hash_idx_mapper[h2];
        if(i1 < i2) {
          bond_order_map[std::make_pair(i1, i2)] = bond_it->GetBondOrder();
        } else {
          bond_order_map[std::make_pair(i2, i1)] = bond_it->GetBondOrder();
        }
      }
    }  
  }

  // super stupid vec intended for sorting... sorry for that
  std::vector<std::pair<std::pair<int, int>, int > > tmp;
  for(auto it = bond_order_map.begin(); it != bond_order_map.end(); ++it) {
    tmp.push_back(std::make_pair(std::make_pair(it->first.first, 
                  it->first.second), it->second));
  }

  std::sort(tmp.begin(), tmp.end());
  for(auto it = tmp.begin(); it != tmp.end(); ++it) {
    bonds.push_back(it->first.first);
    bonds.push_back(it->first.second);
    bond_orders.push_back(it->second);
  }
}

void ResidueDefinition::FromStream(std::istream& stream) {
  Load(stream, name);
  Load(stream, olc);
  Load(stream, chem_type);
  Load(stream, chem_class);
  Load(stream, anames);
  Load(stream, elements);
  LoadIsHetatm(stream, is_hetatm);
  LoadBonds(stream, bonds);
  LoadBondOrders(stream, bond_orders);
}

void ResidueDefinition::ToStream(std::ostream& stream) const{
  Dump(stream, name);
  Dump(stream, olc);
  Dump(stream, chem_type);
  Dump(stream, chem_class);
  Dump(stream, anames);
  Dump(stream, elements);
  DumpIsHetatm(stream, is_hetatm);
  DumpBonds(stream, bonds);
  DumpBondOrders(stream, bond_orders);
}

int ResidueDefinition::GetIdx(const String& aname) const {
  if(idx_mapper.empty()){
    _InitIdxMapper();
  }
  auto it = idx_mapper.find(aname);
  if(it == idx_mapper.end()) {
    return -1;
  } else {
    return it->second;
  }
}

const std::set<int>& ResidueDefinition::GetRotamericAtoms() const {
  if(!rotamer_setup) {
    _InitRotamer();
  }
  return rotameric_atoms;
}

const std::vector<ChiDefinition>& ResidueDefinition::GetChiDefinitions() const {
  if(!rotamer_setup) {
    _InitRotamer();
  }
  return chi_definitions;
}

const std::vector<SidechainAtomRule>&
ResidueDefinition::GetSidechainAtomRules() const {
  if(!rotamer_setup) {
    _InitRotamer();
  }
  return sidechain_atom_rules;
}

int ResidueDefinition::GetNChiAngles() const {
  if(!rotamer_setup) {
    _InitRotamer();
  }
  return chi_definitions.size();
}

void ResidueDefinition::_InitIdxMapper() const {
  idx_mapper.clear();
  for(size_t i = 0; i < anames.size(); ++i) {
    idx_mapper[anames[i]] = i;
  }
}

void ResidueDefinition::_InitRotamer() const {
  rotameric_atoms.clear();
  if(chem_type == 'A') {
    int n_idx = GetIdx("N");
    int ca_idx = GetIdx("CA");
    int cb_idx = GetIdx("CB");
    if(n_idx!=-1 && ca_idx!=-1 && cb_idx!=-1) {
      switch(olc) {
        case 'R': {
          int cg_idx = GetIdx("CG");
          int cd_idx = GetIdx("CD");
          int ne_idx = GetIdx("NE");
          int cz_idx = GetIdx("CZ");
          int nh1_idx = GetIdx("NH1");
          int nh2_idx = GetIdx("NH2");
          if(cg_idx!=-1 && cd_idx!=-1 && ne_idx!=-1 && cz_idx!=-1 &&
             nh1_idx!=-1 && nh2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd_idx);
            _AddChiDefinition(cb_idx, cg_idx, cd_idx, ne_idx);
            _AddChiDefinition(cg_idx, cd_idx, ne_idx, cz_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5475), Real(2.02370926769), 0, Real(0.0));
            _AddAtomRule(cd_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.5384), Real(1.9898498802), 1, Real(0.0));
            _AddAtomRule(ne_idx, cb_idx, cg_idx, cd_idx,
                         Real(1.5034), Real(1.86907309596), 2, Real(0.0));
            _AddAtomRule(cz_idx, cg_idx, cd_idx, ne_idx,
                         Real(1.3401), Real(2.14762764458), 3, Real(0.0));
            _AddAtomRule(nh1_idx, cd_idx, ne_idx, cz_idx,
                         Real(1.3311), Real(2.0605357149), 4, Real(0.0));
            _AddAtomRule(nh2_idx, cd_idx, ne_idx, cz_idx,
                         Real(1.3292), Real(2.13174514839), 4, Real(M_PI));
          }
          break;
        }
        case 'N': {
          int cg_idx = GetIdx("CG");
          int od1_idx = GetIdx("OD1");
          int nd2_idx = GetIdx("ND2");
          if(cg_idx!=-1 && od1_idx!=-1 && nd2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, od1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5319), Real(1.99491133503), 0, Real(0.0));
            _AddAtomRule(od1_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.2323), Real(2.13907553124), 1, Real(0.0));
            _AddAtomRule(nd2_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.3521), Real(2.02719992619), 1, Real(M_PI));
          }
          break;
        }
        case 'D': {
          int cg_idx = GetIdx("CG");
          int od1_idx = GetIdx("OD1");
          int od2_idx = GetIdx("OD2");
          if(cg_idx!=-1 && od1_idx!=-1 && od2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, od1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5218), Real(1.96524073775), 0, Real(0.0));
            _AddAtomRule(od1_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.2565), Real(2.05931398443), 1, Real(0.0));
            _AddAtomRule(od2_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.2541), Real(2.0542525296), 1, Real(M_PI));
          }
          break;
        }
        case 'Q': {
          int cg_idx = GetIdx("CG");
          int cd_idx = GetIdx("CD");
          int oe1_idx = GetIdx("OE1");
          int ne2_idx = GetIdx("NE2");
          if(cg_idx!=-1 && cd_idx!=-1 && oe1_idx!=-1 && ne2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd_idx);
            _AddChiDefinition(cb_idx, cg_idx, cd_idx, oe1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx, 
                         Real(1.5534), Real(2.0162043519), 0, Real(0.0));
            _AddAtomRule(cd_idx, ca_idx, cb_idx, cg_idx, 
                         Real(1.532), Real(1.96349540849), 1, Real(0.0));
            _AddAtomRule(oe1_idx, cb_idx, cg_idx, cd_idx, 
                         Real(1.2294), Real(2.12092410702), 2, Real(0.0));
            _AddAtomRule(ne2_idx, cb_idx, cg_idx, cd_idx, 
                         Real(1.353), Real(2.03924269803), 2, Real(M_PI));
          }
          break;
        }
        case 'E': {
          int cg_idx = GetIdx("CG");
          int cd_idx = GetIdx("CD");
          int oe1_idx = GetIdx("OE1");
          int oe2_idx = GetIdx("OE2");
          if(cg_idx!=-1 && cd_idx!=-1 && oe1_idx!=-1 && oe2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd_idx);
            _AddChiDefinition(cb_idx, cg_idx, cd_idx, oe1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx, 
                         Real(1.5557), Real(2.01917141163), 0, Real(0.0));
            _AddAtomRule(cd_idx, ca_idx, cb_idx, cg_idx, 
                         Real(1.5307), Real(2.01986954333), 1, Real(0.0));
            _AddAtomRule(oe1_idx, cb_idx, cg_idx, cd_idx, 
                         Real(1.259), Real(2.00695410687), 2, Real(0.0));
            _AddAtomRule(oe2_idx, cb_idx, cg_idx, cd_idx, 
                         Real(1.2532), Real(2.09579136579), 2, Real(M_PI));
          }
          break;
        }
        case 'K': {
          int cg_idx = GetIdx("CG");
          int cd_idx = GetIdx("CD");
          int ce_idx = GetIdx("CE");
          int nz_idx = GetIdx("NZ");
          if(cg_idx!=-1 && cd_idx!=-1 && ce_idx!=-1 && nz_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd_idx);
            _AddChiDefinition(cb_idx, cg_idx, cd_idx, ce_idx);
            _AddChiDefinition(cg_idx, cd_idx, ce_idx, nz_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx, 
                         Real(1.5435), Real(2.02039314211), 0, Real(0.0));
            _AddAtomRule(cd_idx, ca_idx, cb_idx, cg_idx, 
                         Real(1.5397), Real(1.97710897666), 1, Real(0.0));
            _AddAtomRule(ce_idx, cb_idx, cg_idx, cd_idx, 
                         Real(1.535), Real(1.96052834877), 2, Real(0.0));
            _AddAtomRule(nz_idx, cg_idx, cd_idx, ce_idx, 
                         Real(1.4604), Real(1.92789069175), 3, Real(0.0));
          }
          break;
        }
        case 'S': {
          int og_idx = GetIdx("OG");
          if(og_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, og_idx);
            _AddAtomRule(og_idx, n_idx, ca_idx, cb_idx,
                         Real(1.4341), Real(1.96262274387), 0, Real(0.0));
          }
          break;
        }
        case 'C': {
          int sg_idx = GetIdx("SG");
          if(sg_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, sg_idx);
            _AddAtomRule(sg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.8359), Real(1.98740641925), 0, Real(0.0));
          }
          break;
        }
        case 'M': {
          int cg_idx = GetIdx("CG");
          int sd_idx = GetIdx("SD");
          int ce_idx = GetIdx("CE");
          if(cg_idx!=-1 && sd_idx!=-1 && ce_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, sd_idx);
            _AddChiDefinition(cb_idx, cg_idx, sd_idx, ce_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.546), Real(2.02318566891), 0, Real(0.0));
            _AddAtomRule(sd_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.8219), Real(1.9247490991), 1, Real(0.0));
            _AddAtomRule(ce_idx, cb_idx, cg_idx, sd_idx,
                         Real(1.8206), Real(1.72682876192), 2, Real(0.0));
          }
          break;
        }
        case 'W': {
          int cg_idx = GetIdx("CG");
          int cd1_idx = GetIdx("CD1");
          int cd2_idx = GetIdx("CD2");
          int ne1_idx = GetIdx("NE1");
          int ce2_idx = GetIdx("CE2");
          int ce3_idx = GetIdx("CE3");
          int cz2_idx = GetIdx("CZ2");
          int cz3_idx = GetIdx("CZ3");
          int ch2_idx = GetIdx("CH2");
          if(cg_idx!=-1 && cd1_idx!=-1 && cd2_idx!=-1 && ne1_idx!=-1 &&
             ce2_idx!=-1 && ce3_idx!=-1 && cz2_idx!=-1 && cz3_idx!=-1 &&
             ch2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5233), Real(2.00957210075), 0, Real(0.0));
            _AddAtomRule(cd1_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.3679), Real(2.25461632773), 1, Real(0.0));
            _AddAtomRule(cd2_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.4407), Real(2.16333560785), 1, Real(M_PI));
            _AddAtomRule(ce2_idx, cd1_idx, cg_idx, cd2_idx,
                         Real(1.4126), Real(1.86139364725), 4, Real(0.0));
            _AddAtomRule(ne1_idx, cg_idx, cd2_idx, ce2_idx,
                         Real(1.3746), Real(1.88268666413), 4, Real(0.0));
            _AddAtomRule(ce3_idx, cd1_idx, cg_idx, cd2_idx,
                         Real(1.4011), Real(2.31325939059), 4, Real(M_PI));
            _AddAtomRule(cz3_idx, ce2_idx, cd2_idx, ce3_idx,
                         Real(1.4017), Real(2.06228104416), 4, Real(0.0));
            _AddAtomRule(ch2_idx, cd2_idx, ce3_idx, cz3_idx,
                         Real(1.4019), Real(2.11132479614), 4, Real(0.0));
            _AddAtomRule(cz2_idx, ce3_idx, cz3_idx, ch2_idx,
                         Real(1.403), Real(2.10957946689), 4, Real(0.0));
          }
          break;
        }
        case 'Y': {
          int cg_idx = GetIdx("CG");
          int cd1_idx = GetIdx("CD1");
          int cd2_idx = GetIdx("CD2");
          int ce1_idx = GetIdx("CE1");
          int ce2_idx = GetIdx("CE2");
          int cz_idx = GetIdx("CZ");
          int oh_idx = GetIdx("OH");
          if(cg_idx!=-1 && cd1_idx!=-1 && cd2_idx!=-1 && ce1_idx!=-1 &&
             ce2_idx!=-1 && cz_idx!=-1 && oh_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5113), Real(1.9711748572), 0, Real(0.0));
            _AddAtomRule(cd1_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.4064), Real(2.10294721573), 1, Real(0.0));
            _AddAtomRule(cd2_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.4068), Real(2.10242361695), 1, Real(M_PI));
            _AddAtomRule(ce1_idx, cd2_idx, cg_idx, cd1_idx,
                         Real(1.4026), Real(2.1013764194), 4, Real(0.0));
            _AddAtomRule(ce2_idx, cd1_idx, cg_idx, cd2_idx,
                         Real(1.4022), Real(2.1041689462), 4, Real(0.0));
            _AddAtomRule(cz_idx, cg_idx, cd1_idx, ce1_idx,
                         Real(1.3978), Real(2.09596589872), 4, Real(0.0));
            _AddAtomRule(oh_idx, cd2_idx, ce2_idx, cz_idx,
                         Real(1.4063), Real(2.09875842552), 4, Real(M_PI));
          }
          break;
        }
        case 'T': {
          int og1_idx = GetIdx("OG1");
          int cg2_idx = GetIdx("CG1");
          if(og1_idx!=-1 && cg2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, og1_idx);
            _AddAtomRule(og1_idx, n_idx, ca_idx, cb_idx,
                         Real(1.4252), Real(1.95756128904), 0, Real(0.0));
            _AddAtomRule(cg2_idx, og1_idx, ca_idx, cb_idx,
                         Real(1.5324), Real(2.02301113599), 4, Real(-2.1665));
          }
          break;
        }

        case 'V': {
          int cg1_idx = GetIdx("CG1");
          int cg2_idx = GetIdx("CG2");
          if(cg1_idx!=-1 && cg2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg1_idx);
            _AddAtomRule(cg1_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5441), Real(1.9891517485), 0, Real(0.0));
            _AddAtomRule(cg2_idx, cg1_idx, ca_idx, cb_idx,
                         Real(1.5414), Real(1.95773582196), 4, Real(2.1640));
          }
          break;
        }
        case 'I': {
          int cg1_idx = GetIdx("CG1");
          int cg2_idx = GetIdx("CG2");
          int cd1_idx = GetIdx("CD1");
          if(cg1_idx!=-1 && cg2_idx!=-1 && cd1_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg1_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg1_idx, cd1_idx);
            _AddAtomRule(cg1_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5498), Real(1.98321762904), 0, Real(0.0));
            _AddAtomRule(cg2_idx, cg1_idx, ca_idx, cb_idx,
                         Real(1.5452), Real(1.9884536168), 4, Real(-2.2696));
            _AddAtomRule(cd1_idx, ca_idx, cb_idx, cg1_idx,
                         Real(1.5381), Real(1.9912461436), 1, Real(0.0));
          }
          break;
        }
        case 'L': {
          int cg_idx = GetIdx("CG");
          int cd1_idx = GetIdx("CD1");
          int cd2_idx = GetIdx("CD2");
          if(cg_idx!=-1 && cd1_idx!=-1 && cd2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5472), Real(2.05006373939), 0, Real(0.0));
            _AddAtomRule(cd1_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.5361), Real(1.9282397576), 1, Real(0.0));
            _AddAtomRule(cd2_idx, cd1_idx, cb_idx, cg_idx,
                         Real(1.536), Real(1.96471713897), 4, Real(2.0944));

          }
          break;
        }
        case 'P': {
          int cg_idx = GetIdx("CG");
          int cd_idx = GetIdx("CD");
          if(cg_idx!=-1 && cd_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5322), Real(1.82194920616), 0, Real(0.0));
            _AddAtomRule(cd_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.5317), Real(1.80135432098), 1, Real(0.0));
          }
          break;
        }
        case 'H': {
          int cg_idx = GetIdx("CG");
          int nd1_idx = GetIdx("ND1");
          int cd2_idx = GetIdx("CD1");
          int ce1_idx = GetIdx("CE2");
          int ne2_idx = GetIdx("NE2");

          if(cg_idx!=-1 && nd1_idx!=-1 && cd2_idx!=-1 && ce1_idx!=-1 &&
             ne2_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, nd1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5109), Real(2.04098802728), 0, Real(0.0));
            _AddAtomRule(nd1_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.3859), Real(2.09736216212), 1, Real(0.0));
            _AddAtomRule(cd2_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.3596), Real(2.26386657276), 1, Real(M_PI));
            _AddAtomRule(ce1_idx, cd2_idx, cg_idx, nd1_idx,
                         Real(1.317), Real(1.8360863731), 4, Real(0.0));
            _AddAtomRule(ne2_idx, nd1_idx, cg_idx, cd2_idx,
                         Real(1.3782), Real(1.84655834861), 4, Real(0.0));
          }
          break;
        }
        case 'F': {
          int cg_idx = GetIdx("CG");
          int cd1_idx = GetIdx("CD1");
          int cd2_idx = GetIdx("CD2");
          int ce1_idx = GetIdx("CE1");
          int ce2_idx = GetIdx("CE2");
          int cz_idx = GetIdx("CZ");
          if(cg_idx!=-1 && cd1_idx!=-1 && cd2_idx!=-1 && ce1_idx!=-1 &&
             ce2_idx!=-1 && cz_idx!=-1) {
            _AddChiDefinition(n_idx, ca_idx, cb_idx, cg_idx);
            _AddChiDefinition(ca_idx, cb_idx, cg_idx, cd1_idx);
            _AddAtomRule(cg_idx, n_idx, ca_idx, cb_idx,
                         Real(1.5109), Real(1.96803326455), 0, Real(0.0));
            _AddAtomRule(cd1_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.4059), Real(2.099980156), 1, Real(0.0));
            _AddAtomRule(cd2_idx, ca_idx, cb_idx, cg_idx,
                         Real(1.4062), Real(2.10765960471), 1, Real(M_PI));
            _AddAtomRule(ce1_idx, cd2_idx, cg_idx, cd1_idx,
                         Real(1.4006), Real(2.10539067668), 4, Real(0.0));
            _AddAtomRule(ce2_idx, cd1_idx, cg_idx, cd2_idx,
                         Real(1.4002), Real(2.10521614376), 4, Real(0.0));
            _AddAtomRule(cz_idx, cg_idx, cd1_idx, ce1_idx,
                         Real(1.4004), Real(2.09317337192), 4, Real(0.0));
          }
          break;
        }
      }
    }
  }

  for(auto it = sidechain_atom_rules.begin();
      it != sidechain_atom_rules.end(); ++it) {
    rotameric_atoms.insert(it->sidechain_atom_idx);
  }

  rotamer_setup = true;
}

void ResidueDefinition::_AddChiDefinition(int idx_one, int idx_two,
                                          int idx_three, int idx_four) const {
  ChiDefinition chi_definition;
  chi_definition.idx_one = idx_one;
  chi_definition.idx_two = idx_two;
  chi_definition.idx_three = idx_three;
  chi_definition.idx_four = idx_four;
  chi_definitions.push_back(chi_definition);

}

void ResidueDefinition::_AddAtomRule(int a_idx, int anch_one_idx,
                                     int anch_two_idx, int anch_three_idx, 
                                     Real bond_length, Real angle,
                                     int dihedral_idx,
                                     Real base_dihedral) const {
  SidechainAtomRule rule;
  rule.sidechain_atom_idx = a_idx;
  rule.anchor_idx[0] = anch_one_idx;
  rule.anchor_idx[1] = anch_two_idx;
  rule.anchor_idx[2] = anch_three_idx;
  rule.bond_length = bond_length;
  rule.angle = angle;
  rule.dihedral_idx = dihedral_idx;
  rule.base_dihedral = base_dihedral;
  sidechain_atom_rules.push_back(rule);
}

BioUnitDefinition::BioUnitDefinition(const ost::io::MMCifInfoBioUnit& bu) {

  au_chains = bu.GetChainList();

  const std::vector<std::pair<int, int> >& bu_ch_intvl =
  bu.GetChainIntervalList();
  for(auto it = bu_ch_intvl.begin(); it != bu_ch_intvl.end(); ++it) {
    chain_intvl.push_back(it->first);      
    chain_intvl.push_back(it->second);      
  }

  const std::vector<std::vector<MMCifInfoTransOpPtr> >& bu_op_list =
  bu.GetOperations();
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

  const std::vector<std::pair<int, int> >& bu_op_intvl =
  bu.GetOperationsIntervalList();
  for(auto it = bu_op_intvl.begin(); it != bu_op_intvl.end(); ++it) {
    op_intvl.push_back(it->first);      
    op_intvl.push_back(it->second);      
  }
}

void BioUnitDefinition::ToStream(std::ostream& stream) const {
  Dump(stream, au_chains);
  Dump(stream, chain_intvl);
  uint32_t size = operations.size();
  stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  for(auto it = operations.begin(); it != operations.end(); ++it) {
    Dump(stream, *it);
  }
  Dump(stream, op_intvl);
}

void BioUnitDefinition::FromStream(std::istream& stream) {
  Load(stream, au_chains);
  Load(stream, chain_intvl);
  uint32_t size = 0;
  stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
  operations.resize(size);
  for(uint i = 0; i < size; ++i) {
    Load(stream, operations[i]);
  }
  Load(stream, op_intvl);
}

ChainData::ChainData(const ost::mol::ChainHandle& chain,
                     const std::vector<ResidueDefinition>& residue_definitions,
                     const std::unordered_map<unsigned long, int>& res_idx_map,
                     const std::vector<std::pair<unsigned long,
                                                 unsigned long> >& 
                     inter_residue_bonds,
                     const std::vector<int>& inter_residue_bond_orders,
                     std::unordered_map<long, int>& atom_idx_mapper) {

  ch_name = chain.GetName();
  chain_type = chain.GetType();

  // process residues
  ost::mol::ResidueHandleList res_list = chain.GetResidueList();
  for(auto res_it = res_list.begin(); res_it != res_list.end(); ++res_it){
    res_def_indices.push_back(res_idx_map.at(res_it->GetHashCode()));
    rnums.push_back(res_it->GetNumber().GetNum());
    insertion_codes.push_back(res_it->GetNumber().GetInsCode());
    sec_structures.push_back(char(res_it->GetSecStructure()));

    const ResidueDefinition& def = residue_definitions[res_def_indices.back()];
    // process atoms of that residue
    for(auto a_it = def.anames.begin(); a_it != def.anames.end(); ++a_it) {
      ost::mol::AtomHandle at = res_it->FindAtom(*a_it);
      atom_idx_mapper[at.GetHashCode()] = positions.size();
      occupancies.push_back(at.GetOccupancy());
      bfactors.push_back(at.GetBFactor());
      positions.push_back(at.GetPos());
    }
  }

  // super stupid vec intended for sorting... sorry for that
  std::vector<std::pair<std::pair<int, int>, int > > tmp;
  for(uint bond_idx = 0; bond_idx < inter_residue_bonds.size(); ++ bond_idx) {
    int at_idx_one = atom_idx_mapper[inter_residue_bonds[bond_idx].first];
    int at_idx_two = atom_idx_mapper[inter_residue_bonds[bond_idx].second];
    if(at_idx_one < at_idx_two) {
      tmp.push_back(std::make_pair(std::make_pair(at_idx_one, at_idx_two), 
                    inter_residue_bond_orders[bond_idx]));
    } else {
      tmp.push_back(std::make_pair(std::make_pair(at_idx_two, at_idx_one), 
                    inter_residue_bond_orders[bond_idx]));

    }
  }

  std::sort(tmp.begin(), tmp.end());
  for(auto it = tmp.begin(); it != tmp.end(); ++it) {
    bonds.push_back(it->first.first);
    bonds.push_back(it->first.second);
    bond_orders.push_back(it->second);
  }
}

void ChainData::ToStream(std::ostream& stream,
                         const std::vector<ResidueDefinition>& res_def,
                         bool lossy, bool avg_bfactors,
                         bool round_bfactors, bool skip_ss,
                         bool infer_aa_pos) const {
  Dump(stream, ch_name);
  if(chain_type > std::numeric_limits<int8_t>::max()) {
    throw ost::Error("ChainType out of bounds");
  }
  int8_t type = chain_type;
  stream.write(reinterpret_cast<char*>(&type), sizeof(int8_t));
  DumpResDefIndices(stream, res_def_indices);
  DumpRnums(stream, rnums);
  DumpInsertionCodes(stream, insertion_codes);
  DumpOccupancies(stream, occupancies);
  if(avg_bfactors) {
    std::vector<Real> tmp;
    int start = 0;
    for(auto it = res_def_indices.begin(); it != res_def_indices.end(); ++it) {
      int len = res_def[*it].anames.size();
      int end = start + len;
      Real avg = 0.0;
      for(int i = start; i < end; ++i) {
        avg += bfactors[i];
      }
      if(len > 0) {
        avg /= len;
      }
      tmp.push_back(avg);
      start += len;
    }
    DumpBFactors(stream, tmp, round_bfactors);
  } else {
    DumpBFactors(stream, bfactors, round_bfactors);
  }

  if(infer_aa_pos) {
    geom::Vec3List positions_to_dump;
    std::vector<Real> chi_angles;
    positions_to_dump.reserve(positions.size());
    int res_start_idx = 0;
    int n_res = res_def_indices.size();
    for(int res_idx = 0; res_idx < n_res; ++res_idx) {
      const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
      std::set<int> skip_indices = def.GetRotamericAtoms();
      int o_idx = -1;
      if(def.chem_type == 'A') {
        // can reconstruct O if there is CA, C,  no OXT, its not the last
        // residue (res_idx < res_def_indices.size()-1) and the next residue is
        // an amino acid too and has N.
        if(def.GetIdx("CA") != -1 && def.GetIdx("C") != -1 &&
           def.GetIdx("OXT") == -1 && res_idx < n_res-1) {
          const ResidueDefinition next_def = res_def[res_def_indices[res_idx+1]];
          if(next_def.chem_type == 'A' && next_def.GetIdx("N") != -1) {
            o_idx = def.GetIdx("O");
          }
        }
      }
      if(o_idx != -1) {
        skip_indices.insert(o_idx);
      }
      int n_atoms = res_def[res_def_indices[res_idx]].anames.size();
      for(int a_idx = 0; a_idx < n_atoms; ++a_idx) {
        // skips atoms in skip_indices
        if(skip_indices.find(a_idx) == skip_indices.end()) {
          positions_to_dump.push_back(positions[res_start_idx+a_idx]);
        }
      }
      const std::vector<ChiDefinition>& chi_vec = def.GetChiDefinitions();
      for(auto it = chi_vec.begin(); it != chi_vec.end(); ++it) {
        Real a = geom::DihedralAngle(positions[res_start_idx + it->idx_one],
                                     positions[res_start_idx + it->idx_two],
                                     positions[res_start_idx + it->idx_three],
                                     positions[res_start_idx + it->idx_four]);
        chi_angles.push_back(a);
      }
      res_start_idx += n_atoms;
    }
    DumpPositions(stream, positions_to_dump, lossy);
    DumpDihedrals(stream, chi_angles);
  } else {
    DumpPositions(stream, positions, lossy);
  }
  DumpBonds(stream, bonds);
  DumpBondOrders(stream, bond_orders);
  if(!skip_ss) {
    DumpSecStructures(stream, sec_structures);
  }
}

void ChainData::FromStream(std::istream& stream,
                           const std::vector<ResidueDefinition>& res_def,
                           int version, bool lossy, bool avg_bfactors,
                           bool round_bfactors, bool skip_ss,
                           bool infer_aa_pos) {
  
  Load(stream, ch_name);
  if(version >= 2) {
    int8_t type;
    stream.read(reinterpret_cast<char*>(&type), sizeof(int8_t));
    chain_type = ost::mol::ChainType(type);
  }
  LoadResDefIndices(stream, res_def_indices);
  LoadRnums(stream, rnums);
  LoadInsertionCodes(stream, insertion_codes);
  LoadOccupancies(stream, occupancies);
  if(avg_bfactors) {
    std::vector<Real> tmp;
    LoadBFactors(stream, tmp, round_bfactors);
    for(size_t i = 0; i < res_def_indices.size(); ++i) {
      int len = res_def[res_def_indices[i]].anames.size();
      Real bfac = tmp[i];
      bfactors.insert(bfactors.end(), len, bfac);
    }
  } else {
    LoadBFactors(stream, bfactors, round_bfactors);
  }
  LoadPositions(stream, positions, lossy);
  if(infer_aa_pos) {
    std::vector<Real> chi_angles;
    LoadDihedrals(stream, chi_angles);

    int n_res = res_def_indices.size();
    int n_at = 0;
    for(auto it = res_def_indices.begin(); it != res_def_indices.end(); ++it) {
      n_at += res_def[*it].anames.size();
    }
    geom::Vec3List full_positions(n_at);

    int pos_idx = 0;
    int full_pos_idx = 0;
    std::vector<bool> infer_rotamer(n_res, false);
    std::vector<bool> infer_oxygen(n_res, false);
    for(int res_idx = 0; res_idx < n_res; ++res_idx) {
      const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
      int n_res_at = def.anames.size();
      std::set<int> inferred_indices = def.GetRotamericAtoms();
      if(!inferred_indices.empty()) {
        infer_rotamer[res_idx] = true;
      }
      if(def.chem_type == 'A') {
        // can reconstruct O if there is CA, C no OXT, its not the last residue
        // (res_idx < res_def.size()-1) and the next residue is an amino acid
        // too and has N.
        int ca_idx = def.GetIdx("CA");
        int c_idx = def.GetIdx("C");
        int o_idx = def.GetIdx("O");
        int oxt_idx = def.GetIdx("OXT");
        if(o_idx != -1 && ca_idx != -1 && c_idx != -1 && oxt_idx == -1 &&
           res_idx < n_res-1) {
          const ResidueDefinition& next = res_def[res_def_indices[res_idx+1]];
          if(next.chem_type == 'A' && next.GetIdx("N") != -1) {
            inferred_indices.insert(o_idx);
            infer_oxygen[res_idx] = true;
          }
        }
      }
      // transfer positions
      for(int i = 0; i < n_res_at; ++i) {
        if(inferred_indices.find(i) == inferred_indices.end()) {
          full_positions[full_pos_idx++] = positions[pos_idx++];
        } else {
          ++full_pos_idx; // skip
        }
      }
    }

    // infer

    int start_idx = 0;
    int chi_idx = 0;
    for(int res_idx = 0; res_idx < n_res; ++res_idx) {
      const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
      int n_res_atoms = def.anames.size();
      if(infer_oxygen[res_idx]) {
        const ResidueDefinition& next_def = res_def[res_def_indices[res_idx+1]];
        int ca_idx = start_idx + def.GetIdx("CA");
        int c_idx = start_idx + def.GetIdx("C");
        int o_idx = start_idx + def.GetIdx("O");
        int n_next_idx = start_idx + n_res_atoms + next_def.GetIdx("N");
        ConstructOPos(full_positions[ca_idx], full_positions[c_idx],
                      full_positions[n_next_idx], full_positions[o_idx]);
      }
      if(infer_rotamer[res_idx]) {
        const std::vector<SidechainAtomRule>& at_rules =
        def.GetSidechainAtomRules();
        std::vector<Real> dihedral_angles;
        for(int i = 0; i < def.GetNChiAngles(); ++i) {
          dihedral_angles.push_back(chi_angles[chi_idx++]);
        }
        for(auto it = at_rules.begin(); it != at_rules.end(); ++it) {
          Real dihedral = it->base_dihedral;
          if(it->dihedral_idx != 4) {
            dihedral += dihedral_angles[it->dihedral_idx];
          }
          ConstructAtomPos(full_positions[start_idx+it->anchor_idx[0]],
                           full_positions[start_idx+it->anchor_idx[1]],
                           full_positions[start_idx+it->anchor_idx[2]],
                           it->bond_length, it->angle, dihedral,
                           full_positions[start_idx+it->sidechain_atom_idx]);
        }
      }
      start_idx += n_res_atoms;
    }
    std::swap(positions, full_positions);
  }
  LoadBonds(stream, bonds);
  LoadBondOrders(stream, bond_orders);
  if(skip_ss) {
    sec_structures.assign(res_def_indices.size(), 'C');
  } else {
    if(version >= 2) {
      LoadSecStructures(stream, sec_structures);
    } else {
      LoadIntVec(stream, sec_structures);
    }
  }
}

DefaultPepLib::DefaultPepLib() {

  /* hardcoded constructor created with:

  from ost import conop
  def ProcessCompound(comp_name, lib, skip_oxt=True, ca_only=False):
      c = lib.FindCompound(comp_name)   
      anames = list()
      idx_mapper = dict()
      element_mapper = dict()
      for a_idx, a in enumerate(c.atom_specs):
          if a.element == "H":
              continue
          if skip_oxt and a.name == "OXT":
              continue
          if ca_only and a.name != "CA":
              continue
          idx_mapper[a_idx] = a.name
          anames.append(a.name)
          element_mapper[a.name] = a.element
      anames.sort()
      bond_data = list()
      for b in c.bond_specs:
          idx_one = b.atom_one
          idx_two = b.atom_two
          if idx_one in idx_mapper and idx_two in idx_mapper:
              aname_one = idx_mapper[idx_one]
              aname_two = idx_mapper[idx_two]
              idx_one = anames.index(aname_one)
              idx_two = anames.index(aname_two)
              if idx_one < idx_two:
                  bond_data.append(((idx_one, idx_two), b.order))
              else:
                  bond_data.append(((idx_two, idx_one), b.order))
      bond_data.sort()
      print(f"  res_def = ResidueDefinition();")
      print(f"  res_def.name = \"{comp_name}\";")
      print(f"  res_def.olc = '{c.GetOneLetterCode()}';")
      print(f"  res_def.chem_type = '{c.chem_type}';")
      print(f"  res_def.chem_class = '{c.chem_class}';")
      for aname in anames:
          print(f"  res_def.anames.push_back(\"{aname}\");")
      for aname in anames:
          print(f"  res_def.elements.push_back(\"{element_mapper[aname]}\");")
      print(f"  res_def.is_hetatm.assign({len(anames)}, false);")
      for b in bond_data:
          print(f"  res_def.bonds.push_back({b[0][0]});")
          print(f"  res_def.bonds.push_back({b[0][1]});")
      for b in bond_data:
          print(f"  res_def.bond_orders.push_back({b[1]});")
      print("  residue_definitions.push_back(res_def);")
      print()
  lib = conop.GetDefaultLib()
  anames = ["ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "LYS", "SER", "CYS",
            "MET", "TRP", "TYR", "THR", "VAL", "ILE", "LEU", "GLY", "PRO",
            "HIS", "PHE"]
  print("  ResidueDefinition res_def;")
  for aname in anames:
    ProcessCompound(aname, lib)
    ProcessCompound(aname, lib, skip_oxt = False)
    ProcessCompound(aname, lib, ca_only=True)
  */
  ResidueDefinition res_def;
  res_def = ResidueDefinition();
  res_def.name = "ALA";
  res_def.olc = 'A';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(5, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(3);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ALA";
  res_def.olc = 'A';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(6, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(3);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ALA";
  res_def.olc = 'A';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ARG";
  res_def.olc = 'R';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("CZ");
  res_def.anames.push_back("N");
  res_def.anames.push_back("NE");
  res_def.anames.push_back("NH1");
  res_def.anames.push_back("NH2");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(11, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(9);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ARG";
  res_def.olc = 'R';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("CZ");
  res_def.anames.push_back("N");
  res_def.anames.push_back("NE");
  res_def.anames.push_back("NH1");
  res_def.anames.push_back("NH2");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(12, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(11);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(9);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ARG";
  res_def.olc = 'R';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ASN";
  res_def.olc = 'N';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("ND2");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OD1");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(8, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ASN";
  res_def.olc = 'N';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("ND2");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OD1");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(9, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ASN";
  res_def.olc = 'N';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ASP";
  res_def.olc = 'D';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OD1");
  res_def.anames.push_back("OD2");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(8, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ASP";
  res_def.olc = 'D';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OD1");
  res_def.anames.push_back("OD2");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(9, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ASP";
  res_def.olc = 'D';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLN";
  res_def.olc = 'Q';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("NE2");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OE1");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(9, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(8);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLN";
  res_def.olc = 'Q';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("NE2");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OE1");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(10, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(8);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLN";
  res_def.olc = 'Q';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLU";
  res_def.olc = 'E';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OE1");
  res_def.anames.push_back("OE2");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(9, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(8);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLU";
  res_def.olc = 'E';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OE1");
  res_def.anames.push_back("OE2");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(10, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(8);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLU";
  res_def.olc = 'E';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "LYS";
  res_def.olc = 'K';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CE");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("NZ");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(9, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "LYS";
  res_def.olc = 'K';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CE");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("NZ");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(10, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "LYS";
  res_def.olc = 'K';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "SER";
  res_def.olc = 'S';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OG");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(6, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "SER";
  res_def.olc = 'S';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OG");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(7, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "SER";
  res_def.olc = 'S';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "CYS";
  res_def.olc = 'C';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("SG");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("S");
  res_def.is_hetatm.assign(6, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "CYS";
  res_def.olc = 'C';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.anames.push_back("SG");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("S");
  res_def.is_hetatm.assign(7, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(6);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "CYS";
  res_def.olc = 'C';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "MET";
  res_def.olc = 'M';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CE");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("SD");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("S");
  res_def.is_hetatm.assign(8, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "MET";
  res_def.olc = 'M';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CE");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.anames.push_back("SD");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("S");
  res_def.is_hetatm.assign(9, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(8);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "MET";
  res_def.olc = 'M';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "TRP";
  res_def.olc = 'W';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CE2");
  res_def.anames.push_back("CE3");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("CH2");
  res_def.anames.push_back("CZ2");
  res_def.anames.push_back("CZ3");
  res_def.anames.push_back("N");
  res_def.anames.push_back("NE1");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(14, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(13);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(11);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(12);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(12);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(10);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "TRP";
  res_def.olc = 'W';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CE2");
  res_def.anames.push_back("CE3");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("CH2");
  res_def.anames.push_back("CZ2");
  res_def.anames.push_back("CZ3");
  res_def.anames.push_back("N");
  res_def.anames.push_back("NE1");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(15, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(13);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(14);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(11);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(12);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(12);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(10);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "TRP";
  res_def.olc = 'W';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "TYR";
  res_def.olc = 'Y';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CE1");
  res_def.anames.push_back("CE2");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("CZ");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OH");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(12, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(11);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "TYR";
  res_def.olc = 'Y';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CE1");
  res_def.anames.push_back("CE2");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("CZ");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OH");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(13, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(12);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(11);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "TYR";
  res_def.olc = 'Y';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "THR";
  res_def.olc = 'T';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CG2");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OG1");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(7, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(6);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "THR";
  res_def.olc = 'T';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CG2");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OG1");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(8, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(6);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "THR";
  res_def.olc = 'T';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "VAL";
  res_def.olc = 'V';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CG1");
  res_def.anames.push_back("CG2");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(7, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "VAL";
  res_def.olc = 'V';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CG1");
  res_def.anames.push_back("CG2");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(8, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "VAL";
  res_def.olc = 'V';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ILE";
  res_def.olc = 'I';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CG1");
  res_def.anames.push_back("CG2");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(8, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ILE";
  res_def.olc = 'I';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CG1");
  res_def.anames.push_back("CG2");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(9, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "ILE";
  res_def.olc = 'I';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "LEU";
  res_def.olc = 'L';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(8, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(5);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "LEU";
  res_def.olc = 'L';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(9, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(5);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "LEU";
  res_def.olc = 'L';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLY";
  res_def.olc = 'G';
  res_def.chem_type = 'A';
  res_def.chem_class = 'P';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(4, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLY";
  res_def.olc = 'G';
  res_def.chem_type = 'A';
  res_def.chem_class = 'P';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(5, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "GLY";
  res_def.olc = 'G';
  res_def.chem_type = 'A';
  res_def.chem_class = 'P';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "PRO";
  res_def.olc = 'P';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(7, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "PRO";
  res_def.olc = 'P';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(8, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "PRO";
  res_def.olc = 'P';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "HIS";
  res_def.olc = 'H';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CE1");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("ND1");
  res_def.anames.push_back("NE2");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(10, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "HIS";
  res_def.olc = 'H';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CE1");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("N");
  res_def.anames.push_back("ND1");
  res_def.anames.push_back("NE2");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(11, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(7);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "HIS";
  res_def.olc = 'H';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "PHE";
  res_def.olc = 'F';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CE1");
  res_def.anames.push_back("CE2");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("CZ");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(11, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(8);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "PHE";
  res_def.olc = 'F';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("C");
  res_def.anames.push_back("CA");
  res_def.anames.push_back("CB");
  res_def.anames.push_back("CD1");
  res_def.anames.push_back("CD2");
  res_def.anames.push_back("CE1");
  res_def.anames.push_back("CE2");
  res_def.anames.push_back("CG");
  res_def.anames.push_back("CZ");
  res_def.anames.push_back("N");
  res_def.anames.push_back("O");
  res_def.anames.push_back("OXT");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("C");
  res_def.elements.push_back("N");
  res_def.elements.push_back("O");
  res_def.elements.push_back("O");
  res_def.is_hetatm.assign(12, false);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(10);
  res_def.bonds.push_back(0);
  res_def.bonds.push_back(11);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(1);
  res_def.bonds.push_back(9);
  res_def.bonds.push_back(2);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(3);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(4);
  res_def.bonds.push_back(7);
  res_def.bonds.push_back(5);
  res_def.bonds.push_back(8);
  res_def.bonds.push_back(6);
  res_def.bonds.push_back(8);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  res_def.bond_orders.push_back(2);
  res_def.bond_orders.push_back(1);
  residue_definitions.push_back(res_def);

  res_def = ResidueDefinition();
  res_def.name = "PHE";
  res_def.olc = 'F';
  res_def.chem_type = 'A';
  res_def.chem_class = 'L';
  res_def.anames.push_back("CA");
  res_def.elements.push_back("C");
  res_def.is_hetatm.assign(1, false);
  residue_definitions.push_back(res_def);
}

OMFPtr OMF::FromEntity(const ost::mol::EntityHandle& ent,
                       uint8_t options) {

  OMFPtr omf(new OMF);
  omf->name_ = ent.GetName();
  omf->options_ = options;
  omf->version_ = OMF_VERSION;

  //////////////////////////////////////////////////////////////////////////////
  // Generate kind of a "mini compound library"... Eeach unique residue gets  //
  // an own entry in the residue_definitions_ vector.                         //
  //////////////////////////////////////////////////////////////////////////////

  std::unordered_map<ResidueDefinition, int, ResidueDefinitionHash> res_def_map;
  std::unordered_map<unsigned long, int> res_idx_map;
  int idx = 0;

  if(omf->OptionSet(DEFAULT_PEPLIB)) {
    omf->residue_definitions_ = DefaultPepLib::Instance().residue_definitions;
    for(auto it = omf->residue_definitions_.begin();
        it != omf->residue_definitions_.end(); ++it) {
      res_def_map[*it] = idx;
      ++idx;
    }
  }

  ost::mol::ResidueHandleList res_list = ent.GetResidueList();
  for(auto it = res_list.begin(); it != res_list.end(); ++it) {
    ResidueDefinition def(*it);
    auto map_it = res_def_map.find(def);
    if(map_it == res_def_map.end()) {
      // we didn't see that one before
      res_def_map[def] = idx;
      res_idx_map[it->GetHashCode()] = idx;
      ++idx;
      omf->residue_definitions_.push_back(def);
    } else {
      res_idx_map[it->GetHashCode()] = map_it->second;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Before processing the single chains, we're dealing with the bonds. Bonds //
  // can only reasonably be extracted from the full EntityHandle or from the  //
  // single atoms. We therefore do some preprocessing here to extract all     //
  // inter-residue bonds of each chain. The intra-residue bonds are dealt     //
  // with in the residue definitions.                                         //
  //////////////////////////////////////////////////////////////////////////////

  ost::mol::ChainHandleList chain_list = ent.GetChainList();
  // get inter-residue bonds for each chain and the corresponding bond orders
  // we use atom hashes here to define the bonds...
  std::vector<std::vector<std::pair<unsigned long, unsigned long> > > 
  inter_residue_bonds(chain_list.size());
  std::vector<std::vector<int> > 
  inter_residue_bond_orders(chain_list.size());

  //TODO we still need a way to deal with interchain bonds!
  std::vector<std::pair<unsigned long, unsigned long> > interchain_bonds;
  std::vector<int> interchain_bond_orders;
  std::vector<std::pair<int, int> > interchain_bond_chain_indices;

  // also keep track of the chain indices
  std::unordered_map<unsigned long, int> chain_idx_map;
  int chain_idx = 0;
  for(auto it = chain_list.begin(); it != chain_list.end(); ++it) {
    chain_idx_map[it->GetHashCode()] = chain_idx++;
  }

  ost::mol::BondHandleList bond_list = ent.GetBondList(); 
  for(auto bond_it = bond_list.begin(); bond_it != bond_list.end(); ++bond_it) {
    const ost::mol::AtomHandle& at_one = bond_it->GetFirst();
    const ost::mol::AtomHandle& at_two = bond_it->GetSecond();
    const ost::mol::ResidueHandle& res_one = at_one.GetResidue();
    const ost::mol::ResidueHandle& res_two = at_two.GetResidue();
    if(res_one.GetChain() == res_two.GetChain()) {
      if(res_one != res_two) {
          if(omf->OptionSet(INFER_PEP_BONDS)) {
              if(res_one.IsPeptideLinking() && res_two.IsPeptideLinking()) {
              String aname_one = at_one.GetName();
              String aname_two = at_two.GetName();
              if((aname_one == "C" && aname_two == "N" &&
                  res_one.GetNext() == res_two) ||
                 (aname_one == "N" && aname_two == "C" &&
                  res_one.GetPrev() == res_two)) {
                Real bond_length = bond_it->GetLength();
                if(bond_length > Real(1.198) && bond_length < Real(1.474)) {
                  // mean bond length (1.336) +- 6 stds (0.023)
                  // this peptide bond can be inferred... skip...
                  continue;
                }
              }
            }
          }
          int idx = chain_idx_map[at_one.GetResidue().GetChain().GetHashCode()];
          inter_residue_bonds[idx].push_back(std::make_pair(at_one.GetHashCode(), 
                                             at_two.GetHashCode()));
          inter_residue_bond_orders[idx].push_back(bond_it->GetBondOrder());
      }
    } else {
      int idx_one = chain_idx_map[at_one.GetResidue().GetChain().GetHashCode()];
      int idx_two = chain_idx_map[at_two.GetResidue().GetChain().GetHashCode()];
      interchain_bonds.push_back(std::make_pair(at_one.GetHashCode(),
                                                at_two.GetHashCode()));
      interchain_bond_orders.push_back(bond_it->GetBondOrder());
      interchain_bond_chain_indices.push_back(std::make_pair(idx_one, idx_two));
    }
  }

  /////////////////////////
  // Fill per chain data //
  /////////////////////////
  // maps for each atom hash code the idx in the respective chain
  std::map<String, std::unordered_map<long, int> > atom_idx_mapper;
  for(uint ch_idx = 0; ch_idx < chain_list.size();  ++ch_idx) {
    String chain_name = chain_list[ch_idx].GetName(); 
    atom_idx_mapper[chain_name] = std::unordered_map<long, int>();
    omf->chain_data_[chain_name] = 
    ChainDataPtr(new ChainData(chain_list[ch_idx], omf->residue_definitions_, 
                 res_idx_map, inter_residue_bonds[ch_idx], 
                 inter_residue_bond_orders[ch_idx],
                 atom_idx_mapper[chain_name]));
  }

  for(uint i = 0; i < interchain_bonds.size(); ++i) {
    int ch_idx_one = interchain_bond_chain_indices[i].first;
    int ch_idx_two = interchain_bond_chain_indices[i].second;
    String chain_name_one = chain_list[ch_idx_one].GetName();
    String chain_name_two = chain_list[ch_idx_two].GetName();
    omf->bond_chain_names_.push_back(chain_name_one);
    omf->bond_chain_names_.push_back(chain_name_two);
    int at_hash_one = interchain_bonds[i].first;
    int at_hash_two = interchain_bonds[i].second;
    omf->bond_atoms_.push_back(atom_idx_mapper[chain_name_one][at_hash_one]);
    omf->bond_atoms_.push_back(atom_idx_mapper[chain_name_two][at_hash_two]);
    omf->bond_orders_.push_back(interchain_bond_orders[i]);
  }

  return omf;
}

OMFPtr OMF::FromMMCIF(const ost::mol::EntityHandle& ent,
                      const MMCifInfo& info,
                      uint8_t options) {

  OMFPtr p = OMF::FromEntity(ent, options);
  const std::vector<MMCifInfoBioUnit>& biounits = info.GetBioUnits();
  for(auto it = biounits.begin(); it != biounits.end(); ++it) {
    p->biounit_definitions_.push_back(BioUnitDefinition(*it));
  }
  return p;
}

OMFPtr OMF::FromFile(const String& fn) {
  std::ifstream in_stream(fn.c_str(), std::ios::binary);
  if (!in_stream) {
    throw ost::Error("Could not open " + fn);
  }
  OMFPtr loaded_omf(new OMF);
  loaded_omf->FromStream(in_stream);
  return loaded_omf;
}

void OMF::ToFile(const String& fn) const {
  std::ofstream out_stream(fn.c_str(), std::ios::binary);
  if (!out_stream) {
    throw ost::Error("Could not open " + fn);
  }
  this->ToStream(out_stream);
}

OMFPtr OMF::FromString(const String& s) {
  std::istringstream in_stream(s);
  OMFPtr loaded_omf(new OMF);
  loaded_omf->FromStream(in_stream);
  return loaded_omf;
}

String OMF::ToString() const {
  std::ostringstream out_stream;
  this->ToStream(out_stream);
  return out_stream.str();
}

ost::mol::EntityHandle OMF::GetAU() const{
  ost::mol::EntityHandle ent = mol::CreateEntity();
  ent.SetName(name_);
  ost::mol::XCSEditor ed = ent.EditXCS(mol::BUFFERED_EDIT);

  for(auto it = chain_data_.begin(); it!=chain_data_.end(); ++it) {
    ost::mol::ChainHandle ch = ed.InsertChain(it->first); 
    this->FillChain(ch, ed, it->second);
  }

  // deal with inter-chain bonds
  // implemented inefficiently but can be considered special case
  for(uint i = 0; i < bond_orders_.size(); ++i) {
    ost::mol::ChainHandle chain_one = ent.FindChain(bond_chain_names_[i*2]);
    ost::mol::ChainHandle chain_two = ent.FindChain(bond_chain_names_[i*2+1]);
    ost::mol::AtomHandle at_one = chain_one.GetAtomList()[bond_atoms_[2*i]];
    ost::mol::AtomHandle at_two = chain_two.GetAtomList()[bond_atoms_[2*i]];
    ed.Connect(at_one, at_two, bond_orders_[i]);
  }
  return ent;
}

ost::mol::EntityHandle OMF::GetAUChain(const String& name) const{

  if(chain_data_.find(name) == chain_data_.end()) {
    throw ost::Error("No chain of name " + name);
  }
  ost::mol::EntityHandle ent = mol::CreateEntity();
  std::stringstream ss;
  ss << name_ << " " << name;
  ent.SetName(ss.str());
  ost::mol::XCSEditor ed = ent.EditXCS(mol::BUFFERED_EDIT);
  ost::mol::ChainHandle ch = ed.InsertChain(name);  
  this->FillChain(ch, ed, chain_data_.at(name));
  return ent;
}

ost::mol::EntityHandle OMF::GetBU(int bu_idx) const{
  if(bu_idx < 0 || bu_idx >= static_cast<int>(biounit_definitions_.size())) {
    throw ost::Error("Invalid biounit idx");
  }

  const BioUnitDefinition& bu = biounit_definitions_[bu_idx];
  ost::mol::EntityHandle ent = mol::CreateEntity();
  std::stringstream ss;
  ss << name_ << " " << bu_idx;
  ent.SetName(ss.str());
  ost::mol::XCSEditor ed = ent.EditXCS(mol::BUFFERED_EDIT);

  std::vector<String> au_chain_names;
  std::vector<geom::Mat4> transforms;

  // The code below is pure magic and heavily inspired by
  // the biounit buildup in modules/io/pymod/__init__.py
  int n_intervals = bu.chain_intvl.size() / 2;
  for(int intvl_idx = 0; intvl_idx < n_intervals; ++intvl_idx) {
    std::vector<geom::Mat4> rts;
    int op_start = bu.op_intvl[2*intvl_idx];
    int op_end = bu.op_intvl[2*intvl_idx+1];
    int n_intv_ops = op_end - op_start;
    if(n_intv_ops) {
      for(auto it = bu.operations[op_start].begin(); 
          it != bu.operations[op_start].end(); ++it) {
        rts.push_back(*it);
      }
      ++op_start;
      while(op_start < op_end) {
        std::vector<geom::Mat4> tmp_rts;
        for(auto i = bu.operations[op_start].begin(); 
            i != bu.operations[op_start].end(); ++i) {
          for(auto j = rts.begin(); j != rts.end(); ++j) {
            tmp_rts.push_back((*j)*(*i));
          }
        }
        rts = tmp_rts;
        ++op_start;
      }
    }
    for(int ch_idx = bu.chain_intvl[2*intvl_idx]; 
        ch_idx < bu.chain_intvl[2*intvl_idx+1]; ++ch_idx) {
      for(auto it = rts.begin(); it != rts.end(); ++it) {
        au_chain_names.push_back(bu.au_chains[ch_idx]);
        transforms.push_back(*it);
      }
    }
  }

  ChainNameGenerator gen;
  for(uint bu_ch_idx = 0; bu_ch_idx < au_chain_names.size(); ++bu_ch_idx) {
    String bu_ch_name = gen.Get();
    ost::mol::ChainHandle added_chain = ed.InsertChain(bu_ch_name);
    this->FillChain(added_chain, ed, chain_data_.at(au_chain_names[bu_ch_idx]),
                    transforms[bu_ch_idx]);
  }

  return ent;
}

void OMF::ToStream(std::ostream& stream) const {

  uint32_t magic_number = 42;
  stream.write(reinterpret_cast<char*>(&magic_number), sizeof(uint32_t));
  // We set it to the current version...
  // If you loaded a structure from a previous version and you dump it again,
  // the version will be updated.
  uint32_t version = version_;
  stream.write(reinterpret_cast<char*>(&version), sizeof(uint32_t));
  stream.write(reinterpret_cast<const char*>(&options_), sizeof(uint8_t));
  DumpName(stream, name_);

  if(OptionSet(DEFAULT_PEPLIB)) {
    // no need to dump the residue definitions from default lib
    auto a = residue_definitions_.begin();
    auto b = residue_definitions_.end();
    int offset = DefaultPepLib::Instance().residue_definitions.size();
    std::vector<ResidueDefinition> tmp(a + offset, b);
    Dump(stream, tmp);
  }
  else {
    Dump(stream, residue_definitions_);
  }

  Dump(stream, biounit_definitions_);
  Dump(stream, chain_data_, residue_definitions_, OptionSet(LOSSY),
       OptionSet(AVG_BFACTORS), OptionSet(ROUND_BFACTORS), OptionSet(SKIP_SS),
       OptionSet(INFER_AA_POS));
  Dump(stream, bond_chain_names_);
  Dump(stream, bond_atoms_);
  Dump(stream, bond_orders_);
}

void OMF::FromStream(std::istream& stream) {

  uint32_t magic_number;
  stream.read(reinterpret_cast<char*>(&magic_number), sizeof(uint32_t));
  if(magic_number != 42) {
    throw ost::Error("Cannot read corrupted OMF stream");
  }

  uint32_t version;
  stream.read(reinterpret_cast<char*>(&version), sizeof(uint32_t));
  if(version != 1 && version != 2) {
    std::stringstream ss;
    ss << "OST version only supports OMF version 1 and 2. Got "<<version;
    throw ost::Error(ss.str());
  }

  version_ = version;

  if(version_ > 1) {
    stream.read(reinterpret_cast<char*>(&options_), sizeof(uint8_t));
    LoadName(stream, name_);
  }

  if(OptionSet(DEFAULT_PEPLIB)) {
    // load residue definitions from default lib and append custom definitions
    std::vector<ResidueDefinition> tmp;
    Load(stream, tmp);
    residue_definitions_ = DefaultPepLib::Instance().residue_definitions;
    residue_definitions_.insert(residue_definitions_.end(),
                                tmp.begin(), tmp.end());
  }
  else {
    Load(stream, residue_definitions_);
  }

  Load(stream, biounit_definitions_);
  Load(stream, chain_data_, residue_definitions_, version_, OptionSet(LOSSY),
       OptionSet(AVG_BFACTORS), OptionSet(ROUND_BFACTORS), OptionSet(SKIP_SS),
       OptionSet(INFER_AA_POS));
  Load(stream, bond_chain_names_);
  Load(stream, bond_atoms_);
  Load(stream, bond_orders_);

  if(!stream.good()) {
    throw ost::Error("Cannot read corrupted OMF stream");
  }
}

void OMF::FillChain(ost::mol::ChainHandle& chain, ost::mol::XCSEditor& ed,
                    const ChainDataPtr data, geom::Mat4 t) const {

  ed.SetChainType(chain, data->chain_type);
  geom::Vec3List* positions = &data->positions;
  geom::Vec3List transformed_positions; // only filled if non-identity transform
  if(t != geom::Mat4()) {
    transformed_positions.resize(positions->size());
    Real a,b,c;
    for (uint i = 0; i < transformed_positions.size(); ++i) {
      const geom::Vec3& p = data->positions[i];
      a = t(0,0)*p[0]+t(0,1)*p[1]+t(0,2)*p[2]+t(0,3);
      b = t(1,0)*p[0]+t(1,1)*p[1]+t(1,2)*p[2]+t(1,3);
      c = t(2,0)*p[0]+t(2,1)*p[1]+t(2,2)*p[2]+t(2,3);
      transformed_positions[i][0] = a;
      transformed_positions[i][1] = b;
      transformed_positions[i][2] = c;
    }
    positions = &transformed_positions; // bend around
  }

  int at_idx = 0;
  for(uint res_idx = 0; res_idx < data->res_def_indices.size(); ++res_idx) {
    const ResidueDefinition& res_def = 
    residue_definitions_[data->res_def_indices[res_idx]];
    ost::mol::ResNum num = ost::mol::ResNum(data->rnums[res_idx], 
                                            data->insertion_codes[res_idx]);
    ost::mol::ResidueHandle res = ed.AppendResidue(chain, res_def.name, num);
    res.SetOneLetterCode(res_def.olc);
    res.SetChemType(ost::mol::ChemType(res_def.chem_type));
    res.SetChemClass(ost::mol::ChemClass(res_def.chem_class));
    res.SetSecStructure(ost::mol::SecStructure(data->sec_structures[res_idx]));

    for(uint i = 0; i < res_def.anames.size(); ++i) {
      ed.InsertAtom(res, res_def.anames[i], (*positions)[at_idx], 
                    res_def.elements[i], data->occupancies[at_idx], 
                    data->bfactors[at_idx], res_def.is_hetatm[i]);
      ++at_idx;
    }

    ost::mol::AtomHandleList added_atoms = res.GetAtomList();
    for(uint bond_idx = 0; bond_idx < res_def.bond_orders.size(); ++bond_idx) {
      ed.Connect(added_atoms[res_def.bonds[2*bond_idx]], 
                 added_atoms[res_def.bonds[2*bond_idx+1]], 
                 res_def.bond_orders[bond_idx]);
    }
  }
  ost::mol::AtomHandleList added_atoms = chain.GetAtomList();
  for(uint bond_idx = 0; bond_idx < data->bond_orders.size(); ++bond_idx) {
    ed.Connect(added_atoms[data->bonds[2*bond_idx]], 
               added_atoms[data->bonds[2*bond_idx+1]], 
               data->bond_orders[bond_idx]);
  }

  if(OptionSet(INFER_PEP_BONDS)) {
    ost::mol::ResidueHandleList res_list = chain.GetResidueList();
    for(size_t i = 1; i < res_list.size(); ++i) {
      if(res_list[i-1].IsPeptideLinking() && res_list[i].IsPeptideLinking()) {
        const ost::mol::AtomHandle& c = res_list[i-1].FindAtom("C");
        const ost::mol::AtomHandle& n = res_list[i].FindAtom("N");
        if(c.IsValid() && n.IsValid()) {
          Real d = geom::Distance(c.GetPos(), n.GetPos());
          if(d > 0.991 && d < 1.681) {
            // mean (1.336) +- 15 stds (0.023)
            // This is an extremely loose threshold but makes sure to also
            // handle inaccuracies that have been introduced with lossy
            // compression
            ed.Connect(c, n);
          }
        }
      }
    }
  }
}

std::vector<String> OMF::GetChainNames() const{
  std::vector<String> chain_names;
  for(auto it = chain_data_.begin(); it != chain_data_.end(); ++it) {
    chain_names.push_back(it->first);
  }
  return chain_names;
}

const geom::Vec3List& OMF::GetPositions(const String& cname) const {
  auto it = chain_data_.find(cname);
  if(it == chain_data_.end()) {
    throw ost::Error("Provided chain name not in OMF structure");
  }
  return it->second->positions;
}

const std::vector<Real>& OMF::GetBFactors(const String& cname) const {
  auto it = chain_data_.find(cname);
  if(it == chain_data_.end()) {
    throw ost::Error("Provided chain name not in OMF structure");
  }
  return it->second->bfactors;
}

std::vector<Real> OMF::GetAvgBFactors(const String& cname) const {
  auto it = chain_data_.find(cname);
  if(it == chain_data_.end()) {
    throw ost::Error("Provided chain name not in OMF structure");
  }
  const std::vector<Real>& bfactors = it->second->bfactors;
  const std::vector<int>& res_def_indices = it->second->res_def_indices;
  std::vector<Real> avg_bfactors;
  avg_bfactors.reserve(it->second->res_def_indices.size());
  int current_atom_idx = 0;
  for(auto i = res_def_indices.begin(); i != res_def_indices.end(); ++i) {
    int size = residue_definitions_[*i].anames.size();
    Real summed_bfac = 0.0;
    for(int j = 0; j < size; ++j) {
      summed_bfac += bfactors[current_atom_idx];
      ++current_atom_idx;
    }
    if(size > 0) {
      summed_bfac /= size;
    }
    avg_bfactors.push_back(summed_bfac);
  }
  return avg_bfactors;
}

String OMF::GetSequence(const String& cname) const {
  auto it = chain_data_.find(cname);
  if(it == chain_data_.end()) {
    throw ost::Error("Provided chain name not in OMF structure");
  }
  const std::vector<int>& indices = it->second->res_def_indices;
  String sequence(indices.size(), 'X');
  for(size_t i = 0; i < indices.size(); ++i) {
    sequence[i] = residue_definitions_[indices[i]].olc;
  }
  return sequence;
}

}} //ns
