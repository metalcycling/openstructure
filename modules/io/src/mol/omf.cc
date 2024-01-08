#include <algorithm>

#include <ost/mol/atom_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/xcs_editor.hh>

#include "omf.hh"

namespace{

  const uint32_t _bit_masks[25] = {
    0,          // 0000 0000 0...
    2147483648, // 1000 0000 0...
    3221225472, // 1100 0000 0...
    3758096384, // 1110 0000 0...
    4026531840, // 1111 0000 0...
    4160749568, // 1111 1000 0...
    4227858432, // 1111 1100 0...
    4261412864, // 1111 1110 0...
    4278190080, // 1111 1111 0...
    4286578688  // 1111 1111 1...
  };

  class BitStorage {
    // stores arbitrary unsigned integers up to 9 bits
    // no range checks performed
  public:
    BitStorage(): buffer_(2, 0), writing_bit_pos_(0), reading_bit_pos_(0) { }

    void Push(uint32_t v, int n_bits) {
      int byte_pos = writing_bit_pos_ / 8;
      v = v << (32 - n_bits - writing_bit_pos_ % 8); // shift the relevant bits to left
      uint32_t tmp = buffer_[byte_pos];
      tmp = tmp << 24; // buffer content at byte_pos occupies first 8 bits
      tmp = tmp | v;
      buffer_[byte_pos] = (tmp >> 24) & 0xFF;
      buffer_[byte_pos+1] = (tmp >> 16) & 0xFF;
      writing_bit_pos_ += n_bits;

      // make sure we have one empty byte in the end
      int new_byte_pos = writing_bit_pos_ / 8;
      if(new_byte_pos > byte_pos) buffer_.push_back(0);
    }

    uint32_t Read(int n_bits) {
      uint32_t reading_mask = _bit_masks[n_bits];
      int right_shift = reading_bit_pos_ % 8;
      reading_mask = reading_mask >> right_shift;
      int byte_pos = reading_bit_pos_ / 8;
      uint32_t tmp = static_cast<uint32_t>(buffer_[byte_pos]) << 24;
      tmp = tmp | (static_cast<uint32_t>(buffer_[byte_pos+1]) << 16);
      reading_bit_pos_ += n_bits;
      return (tmp & reading_mask) >> (32 - n_bits - right_shift);
    }

    void ResetRead() {
      reading_bit_pos_ = 0;
    }

    void Dump(std::ostream& stream) const {
      stream.write(reinterpret_cast<const char*>(&writing_bit_pos_), sizeof(int));
      int n_bytes = std::ceil(static_cast<Real>(writing_bit_pos_) / 8);
      stream.write(reinterpret_cast<const char*>(&buffer_[0]), n_bytes);
    }

    static BitStorage Load(std::istream& stream) {
      BitStorage storage;
      stream.read(reinterpret_cast<char*>(&storage.writing_bit_pos_), sizeof(int));
      int n_bytes = std::ceil(static_cast<Real>(storage.writing_bit_pos_) / 8);
      storage.buffer_.resize(n_bytes + 1, 0);
      stream.read(reinterpret_cast<char*>(&storage.buffer_[0]), n_bytes);
      return storage;
    }

  private:
    std::vector<uint8_t> buffer_;
    int writing_bit_pos_;
    int reading_bit_pos_;
  };

  void ConstructOPos(const geom::Vec3& ca_pos, const geom::Vec3& c_pos,
                     const geom::Vec3& n_next_pos, geom::Vec3& o_pos) {
    geom::Vec3 o_vec = geom::Normalize(geom::Normalize(c_pos - ca_pos) +
                                       geom::Normalize(c_pos - n_next_pos));
    o_pos = c_pos + 1.2339*o_vec;
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

  inline Real PeptideBond() {
    return 1.3310;
  }

  inline Real N_CA_Bond(char olc) {
    switch(olc) {
      case 'G': {
        return 1.4556;
      }
      case 'A': {
        return 1.4618;
      }
      case 'S': {
        return 1.4610;
      }
      case 'C': {
        return 1.4605;
      }
      case 'V': {
        return 1.4614;
      }
      case 'I': {
        return 1.4616;
      }
      case 'L': {
        return 1.4612;
      }
      case 'T': {
        return 1.4606;
      }
      case 'R': {
        return 1.4614;
      }
      case 'K': {
        return 1.4616;
      }
      case 'D': {
        return 1.4624;
      }
      case 'N': {
        return 1.4614;
      }
      case 'E': {
        return 1.4615;
      }
      case 'Q': {
        return 1.4614;
      }
      case 'M': {
        return 1.4618;
      }
      case 'H': {
        return 1.4612;
      }
      case 'P': {
        return 1.4667;
      }
      case 'F': {
        return 1.4608;
      }
      case 'Y': {
        return 1.4609;
      }
      case 'W': {
        return 1.4611;
      }
      default: {
        return 1.46;
      }
    }
  }

  inline Real CA_C_Bond(char olc) {
    switch (olc) {
      case 'G': {
        return 1.5164;
      }
      case 'A': {
        return 1.5255;
      }
      case 'S': {
        return 1.5251;
      }
      case 'C': {
        return 1.5242;
      }
      case 'V': {
        return 1.5258;
      }
      case 'I': {
        return 1.5259;
      }
      case 'L': {
        return 1.5248;
      }
      case 'T': {
        return 1.5254;
      }
      case 'R': {
        return 1.5253;
      }
      case 'K': {
        return 1.5255;
      }
      case 'D': {
        return 1.5268;
      }
      case 'N': {
        return 1.5256;
      }
      case 'E': {
        return 1.5258;
      }
      case 'Q': {
        return 1.5254;
      }
      case 'M': {
        return 1.5249;
      }
      case 'H': {
        return 1.5242;
      }
      case 'P': {
        return 1.5255;
      }
      case 'F': {
        return 1.5243;
      }
      case 'Y': {
        return 1.5242;
      }
      case 'W': {
        return 1.5243;
      }
      default: {
        return 1.52;
      }
    }
  }

  inline Real CA_CB_Bond(char olc) {
    switch(olc) {
      case 'A': {
        return 1.5253;
      }
      case 'S': {
        return 1.5288;
      }
      case 'C': {
        return 1.5294;
      }
      case 'V': {
        return 1.5446;
      }
      case 'I': {
        return 1.5443;
      }
      case 'L': {
        return 1.5310;
      }
      case 'T': {
        return 1.5390;
      }
      case 'R': {
        return 1.5311;
      }
      case 'K': {
        return 1.5312;
      }
      case 'D': {
        return 1.5320;
      }
      case 'N': {
        return 1.5314;
      }
      case 'E': {
        return 1.5316;
      }
      case 'Q': {
        return 1.5308;
      }
      case 'M': {
        return 1.5311;
      }
      case 'H': {
        return 1.5316;
      }
      case 'P': {
        return 1.5332;
      }
      case 'F': {
        return 1.5332;
      }
      case 'Y': {
        return 1.5331;
      }
      case 'W': {
        return 1.5324;
      }
      default: {
        return 1.53;
      }
    }
  }

  inline Real C_CA_CB_Angle(char olc) {
    switch (olc) {
      case 'A': {
        return 1.9255;
      }
      case 'S': {
        return 1.9190;
      }
      case 'C': {
        return 1.9205;
      }
      case 'V': {
        return 1.9299;
      }
      case 'I': {
        return 1.9318;
      }
      case 'L': {
        return 1.9209;
      }
      case 'T': {
        return 1.9170;
      }
      case 'R': {
        return 1.9229;
      }
      case 'K': {
        return 1.9227;
      }
      case 'D': {
        return 1.9247;
      }
      case 'N': {
        return 1.9278;
      }
      case 'E': {
        return 1.9236;
      }
      case 'Q': {
        return 1.9230;
      }
      case 'M': {
        return 1.9215;
      }
      case 'H': {
        return 1.9246;
      }
      case 'P': {
        return 1.9356;
      }
      case 'F': {
        return 1.9236;
      }
      case 'Y': {
        return 1.9225;
      }
      case 'W': {
        return 1.9230;
      }
      default: {
        return 1.92;
      }
    }
  }

  inline Real CA_C_N_Angle(char olc) {
    return 2.0380;
  }

  inline Real C_N_CA_Angle(char olc) {
    return 2.1200;
  }

  inline Real N_CA_C_Angle(char olc) {
    switch(olc) {
      case 'G': {
        return 1.9747;
      }
      case 'A': {
        return 1.9374;
      }
      case 'S': {
        return 1.9401;
      }
      case 'C': {
        return 1.9338;
      }
      case 'V': {
        return 1.9140;
      }
      case 'I': {
        return 1.9139;
      }
      case 'L': {
        return 1.9339;
      }
      case 'T': {
        return 1.9314;
      }
      case 'R': {
        return 1.9359;
      }
      case 'K': {
        return 1.9374;
      }
      case 'D': {
        return 1.9364;
      }
      case 'N': {
        return 1.9434;
      }
      case 'E': {
        return 1.9389;
      }
      case 'Q': {
        return 1.9375;
      }
      case 'M': {
        return 1.9350;
      }
      case 'H': {
        return 1.9374;
      }
      case 'P': {
        return 1.9665;
      }
      case 'F': {
        return 1.9323;
      }
      case 'Y': {
        return 1.9341;
      }
      case 'W': {
        return 1.9348;
      }
      default: {
        return 1.94;
      }
    }
  }


  inline Real N_C_CA_CB_DiAngle(char olc) {
    switch(olc) {
      case 'A': {
        return 2.1423;
      }
      case 'S': {
        return 2.1428;
      }
      case 'C': {
        return 2.1409;
      }
      case 'V': {
        return 2.1533;
      }
      case 'I': {
        return 2.1519;
      }
      case 'L': {
        return 2.1377;
      }
      case 'T': {
        return 2.1496;
      }
      case 'R': {
        return 2.1433;
      }
      case 'K': {
        return 2.1432;
      }
      case 'D': {
        return 2.1449;
      }
      case 'N': {
        return 2.1526;
      }
      case 'E': {
        return 2.1452;
      }
      case 'Q': {
        return 2.1442;
      }
      case 'M': {
        return 2.1419;
      }
      case 'H': {
        return 2.1437;
      }
      case 'P': {
        return 2.0121;
      }
      case 'F': {
        return 2.1426;
      }
      case 'Y': {
        return 2.1423;
      }
      case 'W': {
        return 2.1420;
      }
      default: {
        return 2.14;
      }
    }
  }

  void FillInferredTriPeptideIndices(const ost::io::ResidueDefinition& def_one,
                                     const ost::io::ResidueDefinition& def_two,
                                     const ost::io::ResidueDefinition& def_three,
                                     int res_idx,
                                     std::vector<std::set<int> >& indices) {
    indices[res_idx].insert(def_one.GetIdx("N"));
    indices[res_idx].insert(def_one.GetIdx("C"));
    int cb_idx = def_one.GetIdx("CB");
    if(cb_idx != -1) {
      indices[res_idx].insert(cb_idx);
    }
    indices[res_idx+1].insert(def_two.GetIdx("N"));
    indices[res_idx+1].insert(def_two.GetIdx("C"));
    cb_idx = def_two.GetIdx("CB");
    if(cb_idx != -1) {
      indices[res_idx+1].insert(cb_idx);
    }
    indices[res_idx+2].insert(def_three.GetIdx("N"));
    indices[res_idx+2].insert(def_three.GetIdx("C"));
    cb_idx = def_three.GetIdx("CB");
    if(cb_idx != -1) {
      indices[res_idx+2].insert(cb_idx);
    }
  }

  void FillInferredRotIndices(const ost::io::ResidueDefinition& def,
                              int res_idx,
                              std::vector<std::set<int> >& inferred_indices) {
    const std::vector<ost::io::SidechainAtomRule>& at_rules =
    def.GetSidechainAtomRules();
    for(auto it = at_rules.begin(); it != at_rules.end(); ++it) {
      inferred_indices[res_idx].insert(it->sidechain_atom_idx);
    }
  }

  bool EncodeTriPeptide(const ost::io::ResidueDefinition& def_one,
                        const ost::io::ResidueDefinition& def_two,
                        const ost::io::ResidueDefinition& def_three,
                        Real error_thresh,
                        const std::vector<geom::Vec3>& ref_positions,
                        int res_idx,
                        int res_start_idx,
                        std::vector<std::set<int> >& skip_indices,
                        std::vector<geom::Vec3>& positions,
                        BitStorage& data) {

    // extracts data required to reconstruct positions
    // if max reconstruction error is below specified threshold,
    // the reconstructed positions are directly fed back into
    // positions. skip_indices are updated.

    Real max_error = 0.0;

    int n_one_idx = def_one.GetIdx("N");
    int ca_one_idx = def_one.GetIdx("CA");
    int c_one_idx = def_one.GetIdx("C");
    int n_two_idx = def_two.GetIdx("N");
    int ca_two_idx = def_two.GetIdx("CA");
    int c_two_idx = def_two.GetIdx("C");
    int n_three_idx = def_three.GetIdx("N");
    int ca_three_idx = def_three.GetIdx("CA");
    int c_three_idx = def_three.GetIdx("C");

    if(n_one_idx == -1 || ca_one_idx == -1 || c_one_idx == -1 ||
       n_two_idx == -1 || ca_two_idx == -1 || c_two_idx == -1 ||
       n_three_idx == -1 || ca_three_idx == -1 || c_three_idx == -1) {
      return false;
    }

    // CBeta are optional
    int cb_one_idx = def_one.GetIdx("CB");
    int cb_two_idx = def_two.GetIdx("CB");
    int cb_three_idx = def_three.GetIdx("CB");

    int def_one_size = def_one.anames.size();
    int def_two_size = def_two.anames.size();

    n_one_idx += res_start_idx;
    ca_one_idx += res_start_idx;
    c_one_idx += res_start_idx;
    n_two_idx += (res_start_idx + def_one_size); 
    ca_two_idx += (res_start_idx + def_one_size); 
    c_two_idx += (res_start_idx + def_one_size); 
    n_three_idx += (res_start_idx + def_one_size + def_two_size); 
    ca_three_idx += (res_start_idx + def_one_size + def_two_size); 
    c_three_idx += (res_start_idx + def_one_size + def_two_size);

    if(cb_one_idx != -1) {
      cb_one_idx += res_start_idx;
    }

    if(cb_two_idx != -1) {
      cb_two_idx += (res_start_idx + def_one_size);
    }

    if(cb_three_idx != -1) {
      cb_three_idx += (res_start_idx + def_one_size + def_two_size);
    }

    // derive parameters to reconstruct c_two
    /////////////////////////////////////////
    Real da_c_two = geom::DihedralAngle(positions[ca_one_idx],
                                        positions[ca_three_idx],
                                        positions[ca_two_idx],
                                        ref_positions[c_two_idx]);
    Real a_c_two = geom::Angle(positions[ca_three_idx] - positions[ca_two_idx],
                               ref_positions[c_two_idx] - positions[ca_two_idx]);
    int int_da_c_two = round((da_c_two + M_PI)/(2*M_PI)*255);
    int int_a_c_two = round((a_c_two)/(M_PI)*255);
    geom::Vec3 reconstructed_c_two;
    ConstructAtomPos(positions[ca_one_idx], positions[ca_three_idx],
                     positions[ca_two_idx], CA_C_Bond(def_two.olc),
                     static_cast<Real>(int_a_c_two)/255*M_PI,
                     static_cast<Real>(int_da_c_two)/255*2*M_PI-M_PI,
                     reconstructed_c_two);
    max_error = std::max(max_error, geom::Distance(reconstructed_c_two,
                                                   ref_positions[c_two_idx]));

    // derive parameters to reconstruct n_two
    /////////////////////////////////////////
    Real da_n_two = geom::DihedralAngle(positions[ca_one_idx],
                                        positions[ca_three_idx],
                                        positions[ca_two_idx],
                                        ref_positions[n_two_idx]);
    Real a_n_two = geom::Angle(positions[ca_three_idx] - positions[ca_two_idx],
                               ref_positions[n_two_idx] - positions[ca_two_idx]);
    int int_da_n_two = round((da_n_two + M_PI)/(2*M_PI)*255);
    int int_a_n_two = round((a_n_two)/(M_PI)*255);
    geom::Vec3 reconstructed_n_two;
    ConstructAtomPos(positions[ca_one_idx], positions[ca_three_idx],
                     positions[ca_two_idx], N_CA_Bond(def_two.olc),
                     static_cast<Real>(int_a_n_two)/255*M_PI,
                     static_cast<Real>(int_da_n_two)/255*2*M_PI-M_PI,
                     reconstructed_n_two);
    max_error = std::max(max_error, geom::Distance(reconstructed_n_two,
                                                   ref_positions[n_two_idx]));

    // derive parameters to reconstruct n_three
    ///////////////////////////////////////////
    Real da_n_three = geom::DihedralAngle(reconstructed_n_two,
                                          positions[ca_two_idx],
                                          reconstructed_c_two,
                                          ref_positions[n_three_idx]);
    Real a_n_three = geom::Angle(positions[ca_two_idx] - reconstructed_c_two,
                                 ref_positions[n_three_idx] - reconstructed_c_two);
    int int_da_n_three = round((da_n_three + M_PI)/(2*M_PI)*255);
    // store angle as diff to ideal angle
    Real diff = a_n_three-CA_C_N_Angle(def_two.olc);
    // quantization by 0.5 degrees => 0.0087 in radians
    int int_diff = round(diff/0.0087);
    // make it fit in 4 bits
    int int_a_n_three = std::min(8, std::max(-7, int_diff)) + 7;
    
    geom::Vec3 reconstructed_n_three;
    ConstructAtomPos(reconstructed_n_two,
                     positions[ca_two_idx],
                     reconstructed_c_two, PeptideBond(),
                     CA_C_N_Angle(def_two.olc) + (int_a_n_three-7)*0.0087,
                     static_cast<Real>(int_da_n_three)/255*2*M_PI-M_PI,
                     reconstructed_n_three);
    max_error = std::max(max_error, geom::Distance(reconstructed_n_three,
                                                   ref_positions[n_three_idx]));

    // derive parameters to reconstruct c_three
    ///////////////////////////////////////////
    Real da_c_three = geom::DihedralAngle(reconstructed_c_two,
                                          reconstructed_n_three,
                                          positions[ca_three_idx],
                                          ref_positions[c_three_idx]);
    Real a_c_three = geom::Angle(reconstructed_n_three - positions[ca_three_idx],
                                 ref_positions[c_three_idx] - positions[ca_three_idx]);
    int int_da_c_three = round((da_c_three + M_PI)/(2*M_PI)*255);
    // store angle as diff to ideal angle derived from N_CA_C_Angle function
    diff = a_c_three-N_CA_C_Angle(def_three.olc);
    // quantization by 0.5 degrees => 0.0087 in radians
    int_diff = round(diff/0.0087);
    // make it fit in 4 bits
    int int_a_c_three = std::min(8, std::max(-7, int_diff)) + 7;
    
    geom::Vec3 reconstructed_c_three;
    ConstructAtomPos(reconstructed_c_two,
                     reconstructed_n_three,
                     positions[ca_three_idx], CA_C_Bond(def_three.olc),
                     N_CA_C_Angle(def_three.olc) + (int_a_c_three-7)*0.0087,
                     static_cast<Real>(int_da_c_three)/255*2*M_PI-M_PI,
                     reconstructed_c_three);
    max_error = std::max(max_error, geom::Distance(reconstructed_c_three,
                                                   ref_positions[c_three_idx]));

    // derive parameters to reconstruct c_one
    /////////////////////////////////////////
    Real da_c_one = geom::DihedralAngle(reconstructed_c_two,
                                        positions[ca_two_idx],
                                        reconstructed_n_two,
                                        ref_positions[c_one_idx]);
    Real a_c_one = geom::Angle(positions[ca_two_idx] - reconstructed_n_two,
                               ref_positions[c_one_idx] - reconstructed_n_two);
    int int_da_c_one = round((da_c_one + M_PI)/(2*M_PI)*255);
    // store angle as diff to ideal peptide angle
    diff = a_c_one-C_N_CA_Angle(def_two.olc);
    // quantization by 0.5 degrees => 0.0087 in radians
    int_diff = round(diff/0.0087);
    // make it fit in 4 bits
    int int_a_c_one = std::min(8, std::max(-7, int_diff)) + 7;
    
    geom::Vec3 reconstructed_c_one;
    ConstructAtomPos(reconstructed_c_two,
                     positions[ca_two_idx],
                     reconstructed_n_two, PeptideBond(),
                     C_N_CA_Angle(def_two.olc) + (int_a_c_one-7)*0.0087,
                     static_cast<Real>(int_da_c_one)/255*2*M_PI-M_PI,
                     reconstructed_c_one);
    max_error = std::max(max_error, geom::Distance(reconstructed_c_one,
                                                   ref_positions[c_one_idx]));

    // derive parameters to reconstruct n_one
    /////////////////////////////////////////
    Real da_n_one = geom::DihedralAngle(reconstructed_n_two,
                                        reconstructed_c_one,
                                        positions[ca_one_idx],
                                        ref_positions[n_one_idx]);
    Real a_n_one = geom::Angle(reconstructed_c_one - positions[ca_one_idx],
                               ref_positions[n_one_idx] - positions[ca_one_idx]);
    int int_da_n_one = round((da_n_one + M_PI)/(2*M_PI)*255);
    // store angle as diff to ideal angle derived from N_CA_C_Angle function
    diff = a_n_one-N_CA_C_Angle(def_one.olc);
    // quantization by 0.5 degrees => 0.0087 in radians
    int_diff = round(diff/0.0087);
    // make it fit in 4 bits
    int int_a_n_one = std::min(8, std::max(-7, int_diff)) + 7;
    
    geom::Vec3 reconstructed_n_one;
    ConstructAtomPos(reconstructed_n_two,
                     reconstructed_c_one,
                     positions[ca_one_idx], N_CA_Bond(def_one.olc),
                     N_CA_C_Angle(def_one.olc) + (int_a_n_one-7)*0.0087,
                     static_cast<Real>(int_da_n_one)/255*2*M_PI-M_PI,
                     reconstructed_n_one);
    max_error = std::max(max_error, geom::Distance(reconstructed_n_one,
                                                   ref_positions[n_one_idx]));

    // derive parameters to reconstruct cbetas
    //////////////////////////////////////////
    std::vector<int> cb_data;
    geom::Vec3 reconstructed_cb_one;
    geom::Vec3 reconstructed_cb_two;
    geom::Vec3 reconstructed_cb_three;
    if(cb_one_idx != -1) {
      Real da_cb = geom::DihedralAngle(reconstructed_n_one,
                                       reconstructed_c_one,
                                       positions[ca_one_idx],
                                       ref_positions[cb_one_idx]);
      diff = da_cb - N_C_CA_CB_DiAngle(def_one.olc);
      // quantization by 0.5 degrees => 0.0087 in radians
      int_diff = round(diff/0.0087);
      int int_da_cb = std::min(8, std::max(-7, int_diff)) + 7;

      Real a_cb = geom::Angle(reconstructed_c_one - positions[ca_one_idx],
                              ref_positions[cb_one_idx] - positions[ca_one_idx]);
      diff = a_cb - C_CA_CB_Angle(def_one.olc);
      // quantization by 0.5 degrees => 0.0087 in radians
      int_diff = round(diff/0.0087);
      int int_a_cb = std::min(8, std::max(-7, int_diff)) + 7;

      ConstructAtomPos(reconstructed_n_one,
                       reconstructed_c_one,
                       positions[ca_one_idx], CA_CB_Bond(def_one.olc),
                       C_CA_CB_Angle(def_one.olc) + (static_cast<int>(int_a_cb)-7)*0.0087,
                       N_C_CA_CB_DiAngle(def_one.olc) + (int_da_cb-7) * 0.0087,
                       reconstructed_cb_one);
      max_error = std::max(max_error, geom::Distance(reconstructed_cb_one,
                                                     ref_positions[cb_one_idx]));
      cb_data.push_back(int_a_cb);
      cb_data.push_back(int_da_cb);
    }

    if(cb_two_idx != -1) {
      Real da_cb = geom::DihedralAngle(reconstructed_n_two,
                                       reconstructed_c_two,
                                       positions[ca_two_idx],
                                       ref_positions[cb_two_idx]);
      diff = da_cb - N_C_CA_CB_DiAngle(def_two.olc);
      // quantization by 0.5 degrees => 0.0087 in radians
      int_diff = round(diff/0.0087);
      int int_da_cb = std::min(8, std::max(-7, int_diff)) + 7;

      Real a_cb = geom::Angle(reconstructed_c_two - positions[ca_two_idx],
                              ref_positions[cb_two_idx] - positions[ca_two_idx]);
      diff = a_cb - C_CA_CB_Angle(def_two.olc);
      // quantization by 0.5 degrees => 0.0087 in radians
      int_diff = round(diff/0.0087);
      int int_a_cb = std::min(8, std::max(-7, int_diff)) + 7;

      ConstructAtomPos(reconstructed_n_two,
                       reconstructed_c_two,
                       positions[ca_two_idx], CA_CB_Bond(def_two.olc),
                       C_CA_CB_Angle(def_two.olc) + (int_a_cb-7)*0.0087,
                       N_C_CA_CB_DiAngle(def_two.olc) + (int_da_cb-7) * 0.0087,
                       reconstructed_cb_two);
      max_error = std::max(max_error, geom::Distance(reconstructed_cb_two,
                                                     ref_positions[cb_two_idx]));
      cb_data.push_back(int_a_cb);
      cb_data.push_back(int_da_cb);
    } 

    if(cb_three_idx != -1) {
      Real da_cb = geom::DihedralAngle(reconstructed_n_three,
                                       reconstructed_c_three,
                                       positions[ca_three_idx],
                                       ref_positions[cb_three_idx]);
      diff = da_cb - N_C_CA_CB_DiAngle(def_three.olc);
      // quantization by 0.5 degrees => 0.0087 in radians
      int_diff = round(diff/0.0087);
      int int_da_cb = std::min(8, std::max(-7, int_diff)) + 7;

      Real a_cb = geom::Angle(reconstructed_c_three - positions[ca_three_idx],
                              ref_positions[cb_three_idx] - positions[ca_three_idx]);
      diff = a_cb - C_CA_CB_Angle(def_three.olc);
      // quantization by 0.5 degrees => 0.0087 in radians
      int_diff = round(diff/0.0087);
      int int_a_cb = std::min(8, std::max(-7, int_diff)) + 7;

      ConstructAtomPos(reconstructed_n_three,
                       reconstructed_c_three,
                       positions[ca_three_idx], CA_CB_Bond(def_three.olc),
                       C_CA_CB_Angle(def_three.olc) + (int_a_cb-7)*0.0087,
                       N_C_CA_CB_DiAngle(def_three.olc) + (int_da_cb-7) * 0.0087,
                       reconstructed_cb_three);
      max_error = std::max(max_error, geom::Distance(reconstructed_cb_three,
                                                     ref_positions[cb_three_idx]));
      cb_data.push_back(int_a_cb);
      cb_data.push_back(int_da_cb);
    } 

    if(max_error < error_thresh) {
      positions[n_one_idx] = reconstructed_n_one;
      positions[c_one_idx] = reconstructed_c_one;
      positions[n_two_idx] = reconstructed_n_two;
      positions[c_two_idx] = reconstructed_c_two;
      positions[n_three_idx] = reconstructed_n_three;
      positions[c_three_idx] = reconstructed_c_three;

      if(cb_one_idx != -1) {
        positions[cb_one_idx] = reconstructed_cb_one;
      }

      if(cb_two_idx != -1) {
        positions[cb_two_idx] = reconstructed_cb_two;
      }

      if(cb_three_idx != -1) {
        positions[cb_three_idx] = reconstructed_cb_three;
      }

      FillInferredTriPeptideIndices(def_one, def_two, def_three, res_idx,
                                    skip_indices);

      // push dihedrals to data
      data.Push(int_da_c_two, 8);
      data.Push(int_da_n_two, 8);
      data.Push(int_da_n_three, 8);
      data.Push(int_da_c_three, 8);
      data.Push(int_da_c_one, 8);
      data.Push(int_da_n_one, 8);

      // push angles to data
      data.Push(int_a_c_two, 8);
      data.Push(int_a_n_two, 8);

      // push angle diffs to data
      data.Push(int_a_n_three, 4);
      data.Push(int_a_c_three, 4);
      data.Push(int_a_c_one, 4);
      data.Push(int_a_n_one, 4);

      // push cb data
      for(auto it = cb_data.begin(); it != cb_data.end(); ++it) {
        data.Push(*it, 4);
      }
      return true;
    }
    return false;
  }

  bool EncodePepRotamer(const ost::io::ResidueDefinition& def, Real error_thresh,
                        const std::vector<geom::Vec3>& ref_positions,
                        int res_idx, int res_start_idx,
                        std::vector<std::set<int> >& skip_indices,
                        std::vector<geom::Vec3>& positions,
                        BitStorage& data) {

    const std::vector<ost::io::SidechainAtomRule>& at_rules =
    def.GetSidechainAtomRules();

    int res_n_atoms = def.anames.size();
    std::vector<geom::Vec3> res_ref_positions(ref_positions.begin() + res_start_idx,
                                              ref_positions.begin() + res_start_idx + res_n_atoms);
    std::vector<geom::Vec3> res_positions(positions.begin() + res_start_idx,
                                          positions.begin() + res_start_idx + res_n_atoms);

    // deliberately delay angle computations
    // may lead to tiny corrections for second, third... angle
    // for errors that were introduced for earlier ones
    std::vector<int> comp_dihedrals(def.chi_definitions.size(), 0.0);
    std::vector<bool> dihedral_set(def.chi_definitions.size(), false);
    std::vector<int> angle_diffs;

    for(auto it = at_rules.begin(); it != at_rules.end(); ++it) {
      Real bond = it->bond_length;
      Real angle = it->angle;
      Real dihedral = it->base_dihedral;

      int d_idx = it->dihedral_idx;
      if(d_idx != 5) {
        if(!dihedral_set[d_idx]) {
          const ost::io::ChiDefinition& chi_def = def.chi_definitions[d_idx];
          Real a = geom::DihedralAngle(res_positions[chi_def.idx_one],
                                       res_positions[chi_def.idx_two],
                                       res_positions[chi_def.idx_three],
                                       res_ref_positions[chi_def.idx_four]);
          comp_dihedrals[d_idx] = std::round((a + M_PI)/(2*M_PI)*255);

          dihedral_set[d_idx] = true;
        }
        dihedral += static_cast<Real>(comp_dihedrals[d_idx])/255*2*M_PI-M_PI;
      }

      if(def.critical_sidechain_angles.find(it->sidechain_atom_idx) !=
         def.critical_sidechain_angles.end()) {
        Real a = geom::Angle(res_ref_positions[it->sidechain_atom_idx] -
                             res_positions[it->anchor_idx[2]],
                             res_positions[it->anchor_idx[1]] -
                             res_positions[it->anchor_idx[2]]);
        Real angle_diff = a - angle;
        // quantization by 0.5 degrees => 0.0087 in radians
        int int_diff = std::round(angle_diff/0.0087);
        int_diff = std::min(8, std::max(-7, int_diff)) + 7;
        angle_diffs.push_back(int_diff);

        angle += ((int_diff-7)*0.0087);
      }

      ConstructAtomPos(res_positions[it->anchor_idx[0]],
                       res_positions[it->anchor_idx[1]],
                       res_positions[it->anchor_idx[2]],
                       bond, angle, dihedral,
                       res_positions[it->sidechain_atom_idx]);

      if(geom::Distance(res_positions[it->sidechain_atom_idx],
                        res_ref_positions[it->sidechain_atom_idx]) > error_thresh) {
        return false;
      }
    }

    for(auto it = comp_dihedrals.begin(); it != comp_dihedrals.end(); ++it) {
      data.Push(*it, 8);
    }

    for(auto it = angle_diffs.begin(); it != angle_diffs.end(); ++it) {
      data.Push(*it, 4);
    }

    FillInferredRotIndices(def, res_idx, skip_indices);
    for(auto it = at_rules.begin(); it != at_rules.end(); ++it) {
      positions[res_start_idx + it->sidechain_atom_idx] =
      res_positions[it->sidechain_atom_idx];
    }
    return true;
  }

  void DecodeTriPeptide(const ost::io::ResidueDefinition& def_one,
                        const ost::io::ResidueDefinition& def_two,
                        const ost::io::ResidueDefinition& def_three,
                        int res_start_idx,
                        BitStorage& data,
                        std::vector<geom::Vec3>& positions) {

    int n_one_idx = def_one.GetIdx("N");
    int ca_one_idx = def_one.GetIdx("CA");
    int c_one_idx = def_one.GetIdx("C");
    int n_two_idx = def_two.GetIdx("N");
    int ca_two_idx = def_two.GetIdx("CA");
    int c_two_idx = def_two.GetIdx("C");
    int n_three_idx = def_three.GetIdx("N");
    int ca_three_idx = def_three.GetIdx("CA");
    int c_three_idx = def_three.GetIdx("C");
    int cb_one_idx = def_one.GetIdx("CB");
    int cb_two_idx = def_two.GetIdx("CB");
    int cb_three_idx = def_three.GetIdx("CB");

    int def_one_size = def_one.anames.size();
    int def_two_size = def_two.anames.size();

    n_one_idx += res_start_idx;
    ca_one_idx += res_start_idx;
    c_one_idx += res_start_idx;
    n_two_idx += (res_start_idx + def_one_size);
    ca_two_idx += (res_start_idx + def_one_size);
    c_two_idx += (res_start_idx + def_one_size);
    n_three_idx += (res_start_idx + def_one_size + def_two_size);
    ca_three_idx += (res_start_idx + def_one_size + def_two_size);
    c_three_idx += (res_start_idx + def_one_size + def_two_size);

    if(cb_one_idx != -1) {
      cb_one_idx += res_start_idx;
    }

    if(cb_two_idx != -1) {
      cb_two_idx += (res_start_idx + def_one_size);
    }

    if(cb_three_idx != -1) {
      cb_three_idx += (res_start_idx + def_one_size + def_two_size);
    }

    int int_da_c_two = data.Read(8);
    int int_da_n_two = data.Read(8);
    int int_da_n_three = data.Read(8);
    int int_da_c_three = data.Read(8);
    int int_da_c_one = data.Read(8);
    int int_da_n_one = data.Read(8);

    int int_a_c_two = data.Read(8);
    int int_a_n_two = data.Read(8);

    int int_a_n_three = data.Read(4);
    int int_a_c_three = data.Read(4);

    int int_a_c_one = data.Read(4);
    int int_a_n_one = data.Read(4);

    ConstructAtomPos(positions[ca_one_idx],
                     positions[ca_three_idx],
                     positions[ca_two_idx], CA_C_Bond(def_two.olc),
                     static_cast<Real>(int_a_c_two)/255*M_PI,
                     static_cast<Real>(int_da_c_two)/255*2*M_PI-M_PI,
                     positions[c_two_idx]);

    ConstructAtomPos(positions[ca_one_idx],
                     positions[ca_three_idx],
                     positions[ca_two_idx], N_CA_Bond(def_two.olc),
                     static_cast<Real>(int_a_n_two)/255*M_PI,
                     static_cast<Real>(int_da_n_two)/255*2*M_PI-M_PI,
                     positions[n_two_idx]);

    ConstructAtomPos(positions[n_two_idx],
                     positions[ca_two_idx],
                     positions[c_two_idx], PeptideBond(),
                     CA_C_N_Angle(def_two.olc) + (int_a_n_three-7)*0.0087,
                     static_cast<Real>(int_da_n_three)/255*2*M_PI-M_PI,
                     positions[n_three_idx]);

    ConstructAtomPos(positions[c_two_idx],
                     positions[n_three_idx],
                     positions[ca_three_idx], CA_C_Bond(def_three.olc),
                     N_CA_C_Angle(def_three.olc) + (int_a_c_three-7)*0.0087,
                     static_cast<Real>(int_da_c_three)/255*2*M_PI-M_PI,
                     positions[c_three_idx]);

    ConstructAtomPos(positions[c_two_idx],
                     positions[ca_two_idx],
                     positions[n_two_idx], PeptideBond(),
                     C_N_CA_Angle(def_two.olc) + (int_a_c_one-7)*0.0087,
                     static_cast<Real>(int_da_c_one)/255*2*M_PI-M_PI,
                     positions[c_one_idx]);

    ConstructAtomPos(positions[n_two_idx],
                     positions[c_one_idx],
                     positions[ca_one_idx], N_CA_Bond(def_one.olc),
                     N_CA_C_Angle(def_one.olc) + (int_a_n_one-7)*0.0087,
                     static_cast<Real>(int_da_n_one)/255*2*M_PI-M_PI,
                     positions[n_one_idx]);

    if(cb_one_idx != -1) {
      int int_a_cb = data.Read(4);
      int int_da_cb = data.Read(4);
      ConstructAtomPos(positions[n_one_idx],
                       positions[c_one_idx],
                       positions[ca_one_idx], CA_CB_Bond(def_one.olc),
                       C_CA_CB_Angle(def_one.olc) + (int_a_cb-7)*0.0087,
                       N_C_CA_CB_DiAngle(def_one.olc) + (int_da_cb-7) * 0.0087,
                       positions[cb_one_idx]); 
    }

    if(cb_two_idx != -1) {
      int int_a_cb = data.Read(4);
      int int_da_cb = data.Read(4);
      ConstructAtomPos(positions[n_two_idx],
                       positions[c_two_idx],
                       positions[ca_two_idx], CA_CB_Bond(def_two.olc),
                       C_CA_CB_Angle(def_two.olc) + (int_a_cb-7)*0.0087,
                       N_C_CA_CB_DiAngle(def_two.olc) + (int_da_cb-7) * 0.0087,
                       positions[cb_two_idx]); 
    }

    if(cb_three_idx != -1) {
      int int_a_cb = data.Read(4);
      int int_da_cb = data.Read(4);
      ConstructAtomPos(positions[n_three_idx],
                       positions[c_three_idx],
                       positions[ca_three_idx], CA_CB_Bond(def_three.olc),
                       C_CA_CB_Angle(def_three.olc) + (int_a_cb-7)*0.0087,
                       N_C_CA_CB_DiAngle(def_three.olc) + (int_da_cb-7) * 0.0087,
                       positions[cb_three_idx]); 
    }
  }

  void DecodePepRotamer(const ost::io::ResidueDefinition& def,
                        int res_start_idx, BitStorage& data,
                        std::vector<geom::Vec3>& positions) {
    const std::vector<ost::io::SidechainAtomRule>& at_rules =
    def.GetSidechainAtomRules();

    std::vector<Real> dihedral_angles;
    for(int i = 0; i < def.GetNChiAngles(); ++i) {
      dihedral_angles.push_back(static_cast<Real>(data.Read(8))/255*2*M_PI-M_PI);
    }

    for(auto it = at_rules.begin(); it != at_rules.end(); ++it) {
      Real dihedral = it->base_dihedral;
      if(it->dihedral_idx != 5) {
        dihedral += dihedral_angles[it->dihedral_idx];
      }
      Real angle = it->angle;
      if(def.critical_sidechain_angles.find(it->sidechain_atom_idx) !=
         def.critical_sidechain_angles.end()) {
        int diff = data.Read(4);
        angle += ((diff-7) * 0.0087);
      }
      ConstructAtomPos(positions[res_start_idx+it->anchor_idx[0]],
                       positions[res_start_idx+it->anchor_idx[1]],
                       positions[res_start_idx+it->anchor_idx[2]],
                       it->bond_length, angle, dihedral,
                       positions[res_start_idx+it->sidechain_atom_idx]);
    }
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
            int version, Real max_error, bool avg_bfactors, bool round_bfactors,
            bool skip_ss) {
    uint32_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    map.clear();
    for(uint i = 0; i < size; ++i) {
      ost::io::ChainDataPtr p(new ost::io::ChainData);
      p->FromStream(stream, res_def, version, max_error, avg_bfactors,
                    round_bfactors, skip_ss);
      map[p->ch_name] = p;
    }
  }

  void Dump(std::ostream& stream, 
            const std::map<String, ost::io::ChainDataPtr>& map,
            const std::vector<ost::io::ResidueDefinition>& res_def,
            Real max_error, bool avg_bfactors, bool round_bfactors,
            bool skip_ss) {
    uint32_t size = map.size();
    stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    for(auto it = map.begin(); it != map.end(); ++it) {
        // we don't dump the key (chain name), that's an attribute of the
        // chain itself anyway
      it->second->ToStream(stream, res_def, max_error, avg_bfactors,
                           round_bfactors, skip_ss); 
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

  void Dump(std::ostream& stream,
            const std::vector<bool>& vec) {
    uint32_t size = vec.size();
    uint32_t n_bytes = std::ceil(static_cast<Real>(size)/8);
    std::vector<uint8_t> bit_vector(n_bytes, 0);
    for(size_t i = 0; i < size; ++i) {
      if(vec[i]) bit_vector[i/8] += (1 << (i%8));
    }
    stream.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    stream.write(reinterpret_cast<char*>(&bit_vector[0]),
                                         n_bytes * sizeof(uint8_t));
  }

  void Load(std::istream& stream,
            std::vector<bool>& vec) {
    uint32_t size;
    stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    uint32_t n_bytes = std::ceil(static_cast<Real>(size)/8);
    std::vector<uint8_t> bit_vector(n_bytes);
    stream.read(reinterpret_cast<char*>(&bit_vector[0]),
                                        n_bytes * sizeof(uint8_t));
    vec.resize(size);
    for(uint i = 0; i < size; ++i) {
      vec[i] = static_cast<bool>(bit_vector[i/8] & (1 << (i%8)));
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

  void LoadPositions(std::istream& stream, geom::Vec3List& positions,
                     bool lossy) {

    int8_t n_pos;
    stream.read(reinterpret_cast<char*>(&n_pos), sizeof(int8_t));

    std::vector<Real> x_pos;
    std::vector<Real> y_pos;
    std::vector<Real> z_pos;

    if(n_pos >= 6) {
      int first_x;
      int first_y;
      int first_z;
      stream.read(reinterpret_cast<char*>(&first_x), sizeof(int));
      stream.read(reinterpret_cast<char*>(&first_y), sizeof(int));
      stream.read(reinterpret_cast<char*>(&first_z), sizeof(int));
      std::vector<int> loaded_x;
      std::vector<int> loaded_y;
      std::vector<int> loaded_z;
      Load(stream, loaded_x);
      Load(stream, loaded_y);
      Load(stream, loaded_z);
      std::vector<int> delta_encoded_x;
      std::vector<int> delta_encoded_y;
      std::vector<int> delta_encoded_z;
      delta_encoded_x.push_back(first_x);
      delta_encoded_y.push_back(first_y);
      delta_encoded_z.push_back(first_z);
      delta_encoded_x.insert(delta_encoded_x.end(),
                             loaded_x.begin(), loaded_x.end());
      delta_encoded_y.insert(delta_encoded_y.end(),
                             loaded_y.begin(), loaded_y.end());
      delta_encoded_z.insert(delta_encoded_z.end(),
                             loaded_z.begin(), loaded_z.end());
      std::vector<int> int_x;
      std::vector<int> int_y;
      std::vector<int> int_z;
      DeltaDecoding(delta_encoded_x, int_x);
      DeltaDecoding(delta_encoded_y, int_y);
      DeltaDecoding(delta_encoded_z, int_z);
      if(lossy) {
        IntToRealVec(int_x, x_pos, 0.1);  
        IntToRealVec(int_y, y_pos, 0.1);  
        IntToRealVec(int_z, z_pos, 0.1);  
      } else {
        IntToRealVec(int_x, x_pos, 0.001);  
        IntToRealVec(int_y, y_pos, 0.001);  
        IntToRealVec(int_z, z_pos, 0.001);  
      }
    } else {
      std::vector<int> int_vec;
      Load(stream, int_vec);
      std::vector<Real> real_vec;
      if(lossy) {
        IntToRealVec(int_vec, real_vec, 0.1);
      } else {
        IntToRealVec(int_vec, real_vec, 0.001);
      }
      x_pos.resize(n_pos);
      y_pos.resize(n_pos);
      z_pos.resize(n_pos);
      for(int i = 0; i < n_pos; ++i) {
        x_pos[i] = real_vec[i*3];
        y_pos[i] = real_vec[i*3+1];
        z_pos[i] = real_vec[i*3+2];
      }
    }     

    positions.resize(x_pos.size());
    for(uint i = 0; i < positions.size(); ++i) {
      positions[i] = geom::Vec3(x_pos[i], y_pos[i], z_pos[i]);
    }
  }

  void DumpPositions(std::ostream& stream, const geom::Vec3List& positions,
                     bool lossy) {

    int n_pos = positions.size();

    std::vector<Real> x_pos(n_pos);
    std::vector<Real> y_pos(n_pos);
    std::vector<Real> z_pos(n_pos);
    for(uint i = 0; i < positions.size(); ++i) {
      x_pos[i] = positions[i][0];
      y_pos[i] = positions[i][1];
      z_pos[i] = positions[i][2];
    }

    std::vector<int> int_x;
    std::vector<int> int_y;
    std::vector<int> int_z;

    if(lossy) {
      RealToIntVec(x_pos, int_x, 10);  
      RealToIntVec(y_pos, int_y, 10);  
      RealToIntVec(z_pos, int_z, 10);  
    } else {
      RealToIntVec(x_pos, int_x, 1000);  
      RealToIntVec(y_pos, int_y, 1000);  
      RealToIntVec(z_pos, int_z, 1000);  
    }

    // delta compression is only worth it for a certain amount of
    // positions...
    if(n_pos > 5) {
      // perform delta compression with one quirk: the first values of
      // the delta compressed vectors get dumped explicitely as they
      // may overflow when dumping with small integer sizes and cause
      // excessive integer packing
      int8_t N = 6;
      stream.write(reinterpret_cast<char*>(&N), sizeof(int8_t));
      std::vector<int> x_delta_compressed;
      std::vector<int> y_delta_compressed;
      std::vector<int> z_delta_compressed;
      DeltaEncoding(int_x, x_delta_compressed);
      DeltaEncoding(int_y, y_delta_compressed);
      DeltaEncoding(int_z, z_delta_compressed);
      int first_x = x_delta_compressed[0];
      int first_y = y_delta_compressed[0];
      int first_z = z_delta_compressed[0];
      std::vector<int> x_to_dump(x_delta_compressed.begin() + 1,
                                 x_delta_compressed.end());
      std::vector<int> y_to_dump(y_delta_compressed.begin() + 1,
                                 y_delta_compressed.end());
      std::vector<int> z_to_dump(z_delta_compressed.begin() + 1,
                                 z_delta_compressed.end());
      stream.write(reinterpret_cast<char*>(&first_x), sizeof(int));
      stream.write(reinterpret_cast<char*>(&first_y), sizeof(int));
      stream.write(reinterpret_cast<char*>(&first_z), sizeof(int));
      Dump(stream, x_to_dump);
      Dump(stream, y_to_dump);
      Dump(stream, z_to_dump);
    } else {
      // feed everything in one int vector without any compression and
      // dump it
      int8_t N = n_pos;
      stream.write(reinterpret_cast<char*>(&N), sizeof(int8_t));
      std::vector<int> pos_vec(N * 3);
      for(int i = 0; i < N; ++i) {
        pos_vec[i*3] = int_x[i];
        pos_vec[i*3+1] = int_y[i];
        pos_vec[i*3+2] = int_z[i];
      }
      Dump(stream, pos_vec);
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
  return rotameric_atoms;
}

const std::vector<ChiDefinition>& ResidueDefinition::GetChiDefinitions() const {
  return chi_definitions;
}

const std::vector<SidechainAtomRule>&
ResidueDefinition::GetSidechainAtomRules() const {
  return sidechain_atom_rules;
}

int ResidueDefinition::GetNChiAngles() const {
  return chi_definitions.size();
}

void ResidueDefinition::_InitIdxMapper() const {
  idx_mapper.clear();
  for(size_t i = 0; i < anames.size(); ++i) {
    idx_mapper[anames[i]] = i;
  }
}

void ResidueDefinition::_AddChiDefinition(int idx_one, int idx_two,
                                          int idx_three, int idx_four) {
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
                                     Real base_dihedral) {
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
                         Real max_error, bool avg_bfactors,
                         bool round_bfactors, bool skip_ss) const {
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

  // Lossy means to reduce the accuracy of atom coordinates to one decimal.
  // In terms of eucledian distance, this gives a max error of 0.087. Enable
  // lossy compression if we're above.
  bool lossy = max_error > 0.087;
  // Even when going lower, we might get some lucky shots with internal
  // coordinates. However, at some point it's not worth the overhead...
  bool infer_pos = max_error > 0.05;

  if(infer_pos) {

    int n_res = res_def_indices.size();

    // required info for peptide specific compression
    BitStorage inference_data;
    bool inferred_pep_rotamer = false;
    bool inferred_pep_bb = false;
    bool inferred_pep_o = false;
    std::vector<bool> pep_rotamer_compression;
    std::vector<bool> pep_bb_compression;
    std::vector<bool> pep_o_compression;

    // all indices that can be inferred come in here and won't be dumped
    std::vector<std::set<int> > skip_indices(n_res, std::set<int>());

    // check if we have any peptide residue
    bool pep_present = false;
    for(int res_idx = 0; res_idx < n_res; ++res_idx) {
      const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
      if(def.chem_type == 'A') {
        pep_present = true;
        break;
      }
    }

    if(pep_present) {
      // check for peptide specific compressions
      // same as positions but possibly with reduced accuracy
      std::vector<geom::Vec3> comp_positions = positions;
      if(lossy) {
        for(auto it = comp_positions.begin(); it != comp_positions.end(); ++it) {
          (*it)[0] = 0.1*std::round((*it)[0]*10);
          (*it)[1] = 0.1*std::round((*it)[1]*10);
          (*it)[2] = 0.1*std::round((*it)[2]*10);
        }
      }

      // check tripeptides that can reconstruct with error < 0.5A
      int res_idx = 0;
      int res_start_idx = 0;
      while(res_idx < n_res-2) {

        const ResidueDefinition& res_def_one = res_def[res_def_indices[res_idx]];
        const ResidueDefinition& res_def_two = res_def[res_def_indices[res_idx+1]];
        const ResidueDefinition& res_def_three = res_def[res_def_indices[res_idx+2]];
        if(res_def_one.chem_type == 'A' && res_def_two.chem_type == 'A' &&
           res_def_three.chem_type == 'A') {
          pep_bb_compression.push_back(EncodeTriPeptide(res_def_one,
                                                        res_def_two,
                                                        res_def_three,
                                                        max_error, positions,
                                                        res_idx,
                                                        res_start_idx,
                                                        skip_indices,
                                                        comp_positions,
                                                        inference_data));
          if(pep_bb_compression.back()) {
            res_idx += 3;
            res_start_idx += res_def_one.anames.size();
            res_start_idx += res_def_two.anames.size();
            res_start_idx += res_def_three.anames.size();
            inferred_pep_bb = true;
          } else {
            ++res_idx;
            res_start_idx += res_def_one.anames.size();
          }
        } else {
          // just jump by one residue
          ++res_idx;
          res_start_idx += res_def_one.anames.size();
        }
      }

      // check which residues fulfill 0.5A threshold when applying rotamer
      // compression
      res_start_idx = 0;
      for(res_idx = 0; res_idx < n_res; ++res_idx) {
        const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
        if(def.chem_type == 'A' && !def.GetRotamericAtoms().empty()) {
          pep_rotamer_compression.push_back(EncodePepRotamer(def, max_error,
                                                             positions,
                                                             res_idx,
                                                             res_start_idx,
                                                             skip_indices,
                                                             comp_positions,
                                                             inference_data));
          if(pep_rotamer_compression.back()) {
            inferred_pep_rotamer = true;
          }
          
        }
        res_start_idx += def.anames.size();
      }

      // check which oxygens can be reconstructed within 0.5A
      res_start_idx = 0;
      for(res_idx = 0; res_idx < n_res-1; ++res_idx) {
        const ResidueDefinition& res_def_one = res_def[res_def_indices[res_idx]];
        const ResidueDefinition& res_def_two = res_def[res_def_indices[res_idx+1]];
        int n_atoms = res_def_one.anames.size();
        if(res_def_one.chem_type == 'A' and res_def_two.chem_type == 'A') {
          pep_o_compression.push_back(false);
          int ca_idx = res_def_one.GetIdx("CA");
          int c_idx = res_def_one.GetIdx("C");
          int o_idx = res_def_one.GetIdx("O");
          int oxt_idx = res_def_one.GetIdx("OXT");
          int n_next_idx = res_def_two.GetIdx("N");
          if(ca_idx!=-1 && c_idx!=-1 && o_idx!=-1 && n_next_idx!=-1 &&
             oxt_idx==-1) {
            geom::Vec3 reconstructed_o;
            ConstructOPos(comp_positions[res_start_idx + ca_idx],
                          comp_positions[res_start_idx + c_idx],
                          comp_positions[res_start_idx + n_atoms + n_next_idx],
                          reconstructed_o);
            Real error = geom::Distance(positions[res_start_idx + o_idx],
                              reconstructed_o);
            if(error < max_error) {
              pep_o_compression.back() = true;
              skip_indices[res_idx].insert(o_idx);
              inferred_pep_o = true;
            }
          }
        }
        res_start_idx += n_atoms;
      }
    } // done peptide compression

    int8_t flags = 0;
    if(inferred_pep_rotamer) {
      flags += 1;
    }
    if(inferred_pep_bb) {
      flags += 2;
    }
    if(inferred_pep_o) {
      flags += 4;
    }

    stream.write(reinterpret_cast<char*>(&flags), sizeof(uint8_t));
    if(inferred_pep_rotamer) {
      Dump(stream, pep_rotamer_compression);
    }
    if(inferred_pep_bb) {
      Dump(stream, pep_bb_compression);
    }
    if(inferred_pep_o) {
      Dump(stream, pep_o_compression);
    }
    if(inferred_pep_rotamer || inferred_pep_bb) {
      inference_data.Dump(stream);
    }

    // construct vector containing all positions that cannot be inferred
    geom::Vec3List positions_to_dump;
    int res_start_idx = 0;
    for(int res_idx = 0; res_idx < n_res; ++res_idx) {
      const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
      int res_n_atoms = def.anames.size();

      if(skip_indices[res_idx].empty()) {
        positions_to_dump.insert(positions_to_dump.end(),
                                 positions.begin() + res_start_idx,
                                 positions.begin() + res_start_idx + res_n_atoms);
      } else {
        for(int at_idx = 0; at_idx < res_n_atoms; ++at_idx) {
          if(skip_indices[res_idx].find(at_idx) == skip_indices[res_idx].end()) {
            positions_to_dump.push_back(positions[res_start_idx + at_idx]);
          }
        }
      }
      res_start_idx += res_n_atoms;
    }
    DumpPositions(stream, positions_to_dump, lossy);

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
                           int version, Real max_error, bool avg_bfactors,
                           bool round_bfactors, bool skip_ss) {
  
  Load(stream, ch_name);
  int8_t type;
  stream.read(reinterpret_cast<char*>(&type), sizeof(int8_t));
  chain_type = ost::mol::ChainType(type);

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

  // Lossy means to reduce the accuracy of atom coordinates to one decimal.
  // In terms of eucledian distance, this gives a max error of 0.087. Enable
  // lossy compression if we're above.
  bool lossy = max_error > 0.087;
  // Even when going lower, we might get some lucky shots with internal
  // coordinates. However, at some point it's not worth the overhead...
  bool infer_pos = max_error > 0.05;
  
  if(infer_pos) {

    int n_res = res_def_indices.size();
    int n_at = 0;
    for(auto it = res_def_indices.begin(); it != res_def_indices.end(); ++it) {
      n_at += res_def[*it].anames.size();
    }

    BitStorage inference_data;
    std::vector<bool> pep_bb_compression;
    std::vector<bool> pep_rotamer_compression;
    std::vector<bool> pep_o_compression;

    int8_t flags = 0;
    stream.read(reinterpret_cast<char*>(&flags), sizeof(uint8_t));
    if(flags & 1) {
      Load(stream, pep_rotamer_compression);
    }
    if(flags & 2) {
      Load(stream, pep_bb_compression);
    }
    if(flags & 4) {
      Load(stream, pep_o_compression);
    }
    if(flags & 1 || flags & 2) {
      inference_data = BitStorage::Load(stream);
    }

    // Check which atoms are inferred from the different compression strategies
    std::vector<std::set<int> > inferred_indices(n_res);

    if(!pep_bb_compression.empty()) {
      int res_idx = 0;
      int bb_comp_idx = 0;
      while(res_idx < n_res-2) {
        const ResidueDefinition& res_def_one = res_def[res_def_indices[res_idx]];
        const ResidueDefinition& res_def_two = res_def[res_def_indices[res_idx+1]];
        const ResidueDefinition& res_def_three = res_def[res_def_indices[res_idx+2]];
        if(res_def_one.chem_type == 'A' && res_def_two.chem_type == 'A' &&
           res_def_three.chem_type == 'A') {
          if(pep_bb_compression[bb_comp_idx++]) {
            FillInferredTriPeptideIndices(res_def_one, res_def_two, res_def_three,
                                          res_idx, inferred_indices);
            res_idx += 3;
          } else {
            res_idx += 1;
          }
        }
      }
    }

    if(!pep_rotamer_compression.empty()) {
      int rot_comp_idx = 0;
      for(int res_idx = 0; res_idx < n_res; ++res_idx) {
        const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
        if(def.chem_type == 'A' && !def.GetRotamericAtoms().empty() &&
           pep_rotamer_compression[rot_comp_idx++]) {
          FillInferredRotIndices(def, res_idx, inferred_indices);
        }
      }
    }

    if(!pep_o_compression.empty()) {
      int o_comp_idx = 0;
      for(int res_idx = 0; res_idx < n_res-1; ++res_idx) {
        const ResidueDefinition& res_def_one = res_def[res_def_indices[res_idx]];
        const ResidueDefinition& res_def_two = res_def[res_def_indices[res_idx+1]];
        if(res_def_one.chem_type == 'A' && res_def_two.chem_type == 'A' &&
           pep_o_compression[o_comp_idx++]) {
          inferred_indices[res_idx].insert(res_def_one.GetIdx("O"));
        }
      }
    }

    // fill the positions we have
    LoadPositions(stream, positions, lossy);
    geom::Vec3List full_positions(n_at);

    int pos_idx = 0;
    int full_pos_idx = 0;
    for(int res_idx = 0; res_idx < n_res; ++res_idx) {
      const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
      int n_res_at = def.anames.size();
      if(inferred_indices[res_idx].empty()) {
        for(int i = 0; i < n_res_at; ++i) {
          full_positions[full_pos_idx++] = positions[pos_idx++];
        }
      } else {
        for(int i = 0; i < n_res_at; ++i) {
          if(inferred_indices[res_idx].find(i) == inferred_indices[res_idx].end()) {
            full_positions[full_pos_idx++] = positions[pos_idx++];
          } else {
            ++full_pos_idx; // skip
          }
        }
      }
    }

    // reconstruct the rest
    if(!pep_bb_compression.empty()) {
      int res_idx = 0;
      int bb_comp_idx = 0;
      int res_start_idx = 0;
      while(res_idx < n_res-2) {
        const ResidueDefinition& res_def_one = res_def[res_def_indices[res_idx]];
        const ResidueDefinition& res_def_two = res_def[res_def_indices[res_idx+1]];
        const ResidueDefinition& res_def_three = res_def[res_def_indices[res_idx+2]];
        if(res_def_one.chem_type == 'A' && res_def_two.chem_type == 'A' &&
           res_def_three.chem_type == 'A') {
          if(pep_bb_compression[bb_comp_idx++]) {
            DecodeTriPeptide(res_def_one, res_def_two, res_def_three,
                             res_start_idx, inference_data, full_positions);
            res_idx += 3;
            res_start_idx += res_def_one.anames.size();
            res_start_idx += res_def_two.anames.size();
            res_start_idx += res_def_three.anames.size();
          } else {
            res_idx += 1;
            res_start_idx += res_def_one.anames.size();
          }
        }
      }
    }

    if(!pep_rotamer_compression.empty()) {
      int res_start_idx = 0;
      int rot_comp_idx = 0;
      for(int res_idx = 0; res_idx < n_res; ++res_idx) {
        const ResidueDefinition& def = res_def[res_def_indices[res_idx]];
        if(def.chem_type == 'A' && !def.GetRotamericAtoms().empty()) {
          if(pep_rotamer_compression[rot_comp_idx++]) {
            DecodePepRotamer(def, res_start_idx, inference_data, full_positions);
          }
        }
        res_start_idx += def.anames.size();
      }
    }

    if(!pep_o_compression.empty()) {
      int res_start_idx = 0;
      int o_comp_idx = 0;
      for(int res_idx = 0; res_idx < n_res-1; ++res_idx) {
        const ResidueDefinition& res_def_one = res_def[res_def_indices[res_idx]];
        const ResidueDefinition& res_def_two = res_def[res_def_indices[res_idx+1]];
        int n_atoms = res_def_one.anames.size();
        if(res_def_one.chem_type == 'A' && res_def_two.chem_type == 'A' &&
           pep_o_compression[o_comp_idx++]) {
          int ca_idx = res_def_one.GetIdx("CA");
          int c_idx = res_def_one.GetIdx("C");
          int o_idx = res_def_one.GetIdx("O");
          int n_next_idx = res_def_two.GetIdx("N");
          ConstructOPos(full_positions[res_start_idx + ca_idx],
                        full_positions[res_start_idx + c_idx],
                        full_positions[res_start_idx + n_atoms + n_next_idx],
                        full_positions[res_start_idx + o_idx]);
        }
        res_start_idx += n_atoms;
      }
    }

    std::swap(positions, full_positions);
  } else {
    LoadPositions(stream, positions, lossy);
  }
  LoadBonds(stream, bonds);
  LoadBondOrders(stream, bond_orders);
  if(skip_ss) {
    sec_structures.assign(res_def_indices.size(), 'C');
  } else {
    LoadSecStructures(stream, sec_structures);
  }
}

DefaultPepLib::DefaultPepLib() {

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
  res_def._AddChiDefinition(6, 1, 2, 4);
  res_def._AddChiDefinition(1, 2, 4, 3);
  res_def._AddChiDefinition(2, 4, 3, 7);
  res_def._AddChiDefinition(4, 3, 7, 5);
  res_def._AddChiDefinition(3, 7, 5, 8);
  res_def._AddAtomRule(4, 6, 1, 2, 1.5209, 1.9872, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 4, 1.5220, 1.9507, 1, 0.0); // CD
  res_def._AddAtomRule(7, 2, 4, 3, 1.4604, 1.9486, 2, 0.0); // NE
  res_def._AddAtomRule(5, 4, 3, 7, 1.3304, 2.1771, 3, 0.0); // CZ
  res_def._AddAtomRule(8, 3, 7, 5, 1.3287, 2.1053, 4, 0.0); // NH1
  res_def._AddAtomRule(9, 3, 7, 5, 1.3272, 2.0893, 4, M_PI); // NH2
  res_def.rotameric_atoms.insert({4, 3, 7, 5, 8, 9});
  res_def.critical_sidechain_angles.insert({3, 4, 5, 7});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(4, 1, 2, 3);
  res_def._AddChiDefinition(1, 2, 3, 7);
  res_def._AddAtomRule(3, 4, 1, 2, 1.5154, 1.9661, 0, 0.0); // CG
  res_def._AddAtomRule(7, 1, 2, 3, 1.2329, 2.1098, 1, 0.0); // OD1
  res_def._AddAtomRule(5, 1, 2, 3, 1.3272, 2.0328, 1, M_PI); // ND2
  res_def.rotameric_atoms.insert({3, 7, 5});
  res_def.critical_sidechain_angles.insert({3});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(4, 1, 2, 3);
  res_def._AddChiDefinition(1, 2, 3, 6);
  res_def._AddAtomRule(3, 4, 1, 2, 1.5192, 1.9737, 0, 0.0); // CG
  res_def._AddAtomRule(6, 1, 2, 3, 1.2505, 2.0798, 1, 0.0); // OD1
  res_def._AddAtomRule(7, 1, 2, 3, 1.2508, 2.0593, 1, M_PI); // OD2
  res_def.rotameric_atoms.insert({3, 6, 7});
  res_def.critical_sidechain_angles.insert({3});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(5, 1, 2, 4);
  res_def._AddChiDefinition(1, 2, 4, 3);
  res_def._AddChiDefinition(2, 4, 3, 8);
  res_def._AddAtomRule(4, 5, 1, 2, 1.52, 1.9853, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 4, 1.52, 1.9684, 1, 0.0); // CD
  res_def._AddAtomRule(8, 2, 4, 3, 1.24, 2.1094, 2, 0.0); // OE1
  res_def._AddAtomRule(6, 2, 4, 3, 1.33, 2.0333, 2, M_PI); // NE2
  res_def.rotameric_atoms.insert({4, 3, 8, 6});
  res_def.critical_sidechain_angles.insert({3, 4});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(5, 1, 2, 4);
  res_def._AddChiDefinition(1, 2, 4, 3);
  res_def._AddChiDefinition(2, 4, 3, 7);
  res_def._AddAtomRule(4, 5, 1, 2, 1.5215, 1.9869, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 4, 1.5212, 1.9786, 1, 0.0); // CD
  res_def._AddAtomRule(7, 2, 4, 3, 1.2522, 2.0762, 2, 0.0); // OE1
  res_def._AddAtomRule(8, 2, 4, 3, 1.2516, 2.0610, 2, M_PI); // OE2
  res_def.rotameric_atoms.insert({4, 3, 7, 8});
  res_def.critical_sidechain_angles.insert({3, 4});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(6, 1, 2, 5);
  res_def._AddChiDefinition(1, 2, 5, 3);
  res_def._AddChiDefinition(2, 5, 3, 4);
  res_def._AddChiDefinition(5, 3, 4, 7);
  res_def._AddAtomRule(5, 6, 1, 2, 1.5217, 1.9903, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 5, 1.5230, 1.9488, 1, 0.0); // CD
  res_def._AddAtomRule(4, 2, 5, 3, 1.5215, 1.9493, 2, 0.0); // CE
  res_def._AddAtomRule(7, 5, 3, 4, 1.4922, 1.9498, 3, 0.0); // NZ
  res_def.rotameric_atoms.insert({5, 3, 4, 7});
  res_def.critical_sidechain_angles.insert({3, 4, 5});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(3, 1, 2, 5);
  res_def._AddAtomRule(5, 3, 1, 2, 1.4171, 1.9335, 0, 0.0); // OG
  res_def.rotameric_atoms.insert(5);
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
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
  res_def._AddChiDefinition(3, 1, 2, 5);
  res_def._AddAtomRule(5, 3, 1, 2, 1.8072, 1.9860, 0, 0.0); // SG
  res_def.rotameric_atoms.insert(5);
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
  res_def._AddChiDefinition(3, 1, 2, 6);
  res_def._AddAtomRule(6, 3, 1, 2, 1.808, 1.9865, 0, 0.0); // SG
  res_def.rotameric_atoms.insert(6);
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
  res_def._AddChiDefinition(5, 1, 2, 4);
  res_def._AddChiDefinition(1, 2, 4, 7);
  res_def._AddChiDefinition(2, 4, 7, 3);
  res_def._AddAtomRule(4, 5, 1, 2, 1.5199, 1.9856, 0, 0.0); // CG
  res_def._AddAtomRule(7, 1, 2, 4, 1.8063, 1.9668, 1, 0.0); // SD
  res_def._AddAtomRule(3, 2, 4, 7, 1.7868, 1.7561, 2, 0.0); // CE
  res_def.rotameric_atoms.insert({4, 7, 3});
  res_def.critical_sidechain_angles.insert({4, 7});
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
  res_def._AddChiDefinition(5, 1, 2, 4);
  res_def._AddChiDefinition(1, 2, 4, 8);
  res_def._AddChiDefinition(2, 4, 8, 3);
  res_def._AddAtomRule(4, 5, 1, 2, 1.5199, 1.9856, 0, 0.0); // CG
  res_def._AddAtomRule(8, 1, 2, 4, 1.8063, 1.9668, 1, 0.0); // SD
  res_def._AddAtomRule(3, 2, 4, 8, 1.7868, 1.7561, 2, 0.0); // CE
  res_def.rotameric_atoms.insert({4, 8, 3});
  res_def.critical_sidechain_angles.insert({4, 8});
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
  res_def._AddChiDefinition(11, 1, 2, 7);
  res_def._AddChiDefinition(1, 2, 7, 3);
  res_def._AddAtomRule(7, 11, 1, 2, 1.4986, 1.9896, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 7,  1.3674, 2.2179, 1, 0.0); // CD1
  res_def._AddAtomRule(4, 1, 2, 7,  1.4330, 2.2097, 1, M_PI); // CD2
  res_def._AddAtomRule(5, 3, 7, 4,  1.4127, 1.8710, 5, 0.0); // CE2
  res_def._AddAtomRule(12, 2, 7, 3, 1.3749, 1.9219, 5, M_PI); // NE1 
  res_def._AddAtomRule(6, 3, 7, 4,  1.4001, 2.3370, 5, M_PI); // CE3
  res_def._AddAtomRule(10, 5, 4, 6, 1.3882, 2.0715, 5, 0.0); // CZ3
  res_def._AddAtomRule(8, 4, 6, 10, 1.4025, 2.1130, 5, 0.0); // CH2
  res_def._AddAtomRule(9, 6, 10, 8, 1.3714, 2.1213, 5, 0.0); // CZ2
  res_def.rotameric_atoms.insert({7, 3, 4, 5, 12, 6, 10, 8, 9});
  res_def.critical_sidechain_angles.insert({3, 4, 6, 7, 8, 10});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(9, 1, 2, 7);
  res_def._AddChiDefinition(1, 2, 7, 3);
  res_def._AddAtomRule(7, 9, 1, 2, 1.5104, 1.9853, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 7, 1.3910, 2.1111, 1, 0.0); // CD1
  res_def._AddAtomRule(4, 1, 2, 7, 1.3903, 2.1093, 1, M_PI); // CD2
  res_def._AddAtomRule(5, 4, 7, 3, 1.3888, 2.1144, 5, 0.0); // CE1
  res_def._AddAtomRule(6, 3, 7, 4, 1.3885, 2.1147, 5, 0.0); // CE2
  res_def._AddAtomRule(8, 7, 3, 5, 1.3814, 2.0866, 5, 0.0); // CZ
  res_def._AddAtomRule(11, 3, 5, 8, 1.3771, 2.0909, 5, M_PI); // OH
  res_def.rotameric_atoms.insert({7, 3, 4, 5, 6, 8, 11});
  res_def.critical_sidechain_angles.insert({3, 4, 5, 7, 8});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(4, 1, 2, 6);
  res_def._AddAtomRule(6, 4, 1, 2, 1.4323, 1.9059, 0, 0.0); // OG1
  res_def._AddAtomRule(3, 6, 1, 2, 1.5239, 1.9412, 5, -2.1068); // CG2
  res_def.rotameric_atoms.insert({6, 3});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
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
  res_def._AddChiDefinition(5, 1, 2, 3);
  res_def._AddAtomRule(3, 5, 1, 2, 1.5262, 1.9338, 0, 0.0); // CG1
  res_def._AddAtomRule(4, 3, 1, 2, 1.5257, 1.9294, 5, 2.1478); // CG2
  res_def.rotameric_atoms.insert({3, 4});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
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
  res_def._AddChiDefinition(6, 1, 2, 4);
  res_def._AddChiDefinition(1, 2, 4, 3);
  res_def._AddAtomRule(4, 6, 1, 2, 1.5331, 1.9286, 0, 0.0); // CG1
  res_def._AddAtomRule(5, 4, 1, 2, 1.5300, 1.9320, 5, -2.1534); // CG2
  res_def._AddAtomRule(3, 1, 2, 4, 1.5194, 1.9887, 1, 0.0); // CD1
  res_def.rotameric_atoms.insert({4, 5, 3});
  res_def.critical_sidechain_angles.insert({4});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(6, 1, 2, 5);
  res_def._AddChiDefinition(1, 2, 5, 3);
  res_def._AddAtomRule(5, 6, 1, 2, 1.5298, 2.0276, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 5, 1.5236, 1.9276, 1, 0.0); // CD1
  res_def._AddAtomRule(4, 3, 2, 5, 1.5242, 1.9316, 5, 2.1459); // CD2
  res_def.rotameric_atoms.insert({5, 3, 4});
  res_def.critical_sidechain_angles.insert({5});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(5, 1, 2, 4);
  res_def._AddChiDefinition(1, 2, 4, 3);
  res_def._AddAtomRule(4, 5, 1, 2, 1.4955, 1.8224, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 4, 1.5063, 1.8391, 1, 0.0); // CD
  res_def.rotameric_atoms.insert({4, 3});
  res_def.critical_sidechain_angles.insert({4});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
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
  res_def._AddChiDefinition(6, 1, 2, 5);
  res_def._AddChiDefinition(1, 2, 5, 7);
  res_def._AddAtomRule(5, 6, 1, 2, 1.4953, 1.9832, 0, 0.0); // CG
  res_def._AddAtomRule(7, 1, 2, 5, 1.3783, 2.1399, 1, 0.0); // ND1
  res_def._AddAtomRule(3, 1, 2, 5, 1.3551, 2.2869, 1, M_PI); // CD2
  res_def._AddAtomRule(4, 3, 5, 7, 1.3234, 1.9046, 5, 0.0); // CE1
  res_def._AddAtomRule(8, 7, 5, 3, 1.3734, 1.8712, 5, 0.0); // NE2
  res_def.rotameric_atoms.insert({5, 7, 3, 4, 8});
  res_def.critical_sidechain_angles.insert({3, 5, 7});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
  res_def._AddChiDefinition(9, 1, 2, 7);
  res_def._AddChiDefinition(1, 2, 7, 3);
  res_def._AddAtomRule(7, 9, 1, 2, 1.5043, 1.9866, 0, 0.0); // CG
  res_def._AddAtomRule(3, 1, 2, 7, 1.3883, 2.1071, 1, 0.0); // CD1
  res_def._AddAtomRule(4, 1, 2, 7, 1.3879, 2.1039, 1, M_PI); // CD2
  res_def._AddAtomRule(5, 4, 7, 3, 1.3906, 2.1079, 5, 0.0); // CE1
  res_def._AddAtomRule(6, 3, 7, 4, 1.3902, 2.1082, 5, 0.0); // CE2
  res_def._AddAtomRule(8, 7, 3, 5, 1.3832, 2.0928, 5, 0.0); // CZ
  res_def.rotameric_atoms.insert({7, 3, 4, 5, 6, 8});
  res_def.critical_sidechain_angles.insert({3, 4, 5, 7});
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
  // same rotamer information as previous one
  res_def.chi_definitions = residue_definitions.back().chi_definitions;
  res_def.sidechain_atom_rules = residue_definitions.back().sidechain_atom_rules;
  res_def.rotameric_atoms = residue_definitions.back().rotameric_atoms;
  res_def.critical_sidechain_angles = residue_definitions.back().critical_sidechain_angles;
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
                       Real max_error,
                       uint8_t options) {

  if(max_error < 0.0 || max_error > 10.0) {
    throw ost::Error("max_error must be in [0.0, 10.0]");
  }

  OMFPtr omf(new OMF);
  omf->name_ = ent.GetName();
  omf->max_error_ = 1000 * max_error;
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
    this->FillChain(it->second, ed, ch);
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
  this->FillChain(chain_data_.at(name), ed, ch);
  return ent;
}

void OMF::ToStream(std::ostream& stream) const {

  uint16_t magic_number = 42;
  stream.write(reinterpret_cast<char*>(&magic_number), sizeof(uint16_t));
  // We set it to the current version...
  // If you loaded a structure from a previous version and you dump it again,
  // the version will be updated.
  uint8_t version = version_;
  stream.write(reinterpret_cast<char*>(&version), sizeof(uint8_t));
  stream.write(reinterpret_cast<const char*>(&max_error_), sizeof(uint16_t));
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

  Dump(stream, chain_data_, residue_definitions_, this->GetMaxError(),
       OptionSet(AVG_BFACTORS), OptionSet(ROUND_BFACTORS), OptionSet(SKIP_SS));
  Dump(stream, bond_chain_names_);
  Dump(stream, bond_atoms_);
  Dump(stream, bond_orders_);
}

void OMF::FromStream(std::istream& stream) {

  uint16_t magic_number;
  stream.read(reinterpret_cast<char*>(&magic_number), sizeof(uint16_t));
  if(magic_number != 42) {
    throw ost::Error("Cannot read corrupted OMF stream");
  }

  uint8_t version;
  stream.read(reinterpret_cast<char*>(&version), sizeof(uint8_t));
  if(version < 3) {
    std::stringstream ss;
    ss << "Old OMF versions are deprecated. Can only load versions >= 3, ";
    ss << "got "<< static_cast<int>(version);
    throw ost::Error(ss.str());
  }

  version_ = version;

  stream.read(reinterpret_cast<char*>(&max_error_), sizeof(uint16_t));
  stream.read(reinterpret_cast<char*>(&options_), sizeof(uint8_t));
  LoadName(stream, name_);

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

  Load(stream, chain_data_, residue_definitions_, version_, this->GetMaxError(),
       OptionSet(AVG_BFACTORS), OptionSet(ROUND_BFACTORS), OptionSet(SKIP_SS));
  Load(stream, bond_chain_names_);
  Load(stream, bond_atoms_);
  Load(stream, bond_orders_);

  if(!stream.good()) {
    throw ost::Error("Cannot read corrupted OMF stream");
  }
}

void OMF::FillChain(const ChainDataPtr data, ost::mol::XCSEditor& ed,
                    ost::mol::ChainHandle& chain) const {

  ed.SetChainType(chain, data->chain_type);
  const geom::Vec3List& positions = data->positions;

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
      ed.InsertAtom(res, res_def.anames[i], positions[at_idx], 
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
