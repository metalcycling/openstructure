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

/*
  Author: Gerardo Tauriello, Gabriel Studer
 */

#ifndef OST_SEQ_PROFILE_HANDLE_HH
#define OST_SEQ_PROFILE_HANDLE_HH

#include <ost/base.hh>
#include <ost/stdint.hh>
#include <ost/message.hh>
#include <ost/seq/module_config.hh>

#include <string.h> // for memcpy, etc
#include <vector>
#include <map>
#include <fstream>
#include <boost/shared_ptr.hpp>

namespace ost { namespace seq { 

class ProfileHandle;
class ProfileColumn;
class HMMData;
class ProfileDB;
typedef boost::shared_ptr<ProfileHandle> ProfileHandlePtr;
typedef boost::shared_ptr<HMMData> HMMDataPtr;
typedef std::vector<ProfileHandlePtr> ProfileHandleList;
typedef boost::shared_ptr<ProfileDB> ProfileDBPtr;
typedef std::vector<ProfileColumn> ProfileColumnList;

typedef enum {
  HMM_M2M = 0, HMM_M2I = 1, HMM_M2D = 2,
  HMM_I2M = 3, HMM_I2I = 4,
  HMM_D2M = 5, HMM_D2D = 6
} HMMTransition;

class DLLEXPORT_OST_SEQ HMMData {
public:
  HMMData() {
    memset(trans_, 0, 7*sizeof(Real));
    trans_[HMM_M2M] = 1.0;
    trans_[HMM_I2M] = 1.0;
    trans_[HMM_D2M] = 1.0;
    neff_ = 1.0;
    neff_i_ = 1.0;
    neff_d_ = 1.0;
  }

  HMMData(const HMMData& rhs) {
    memcpy(trans_, rhs.trans_, 7*sizeof(Real));
    neff_ = rhs.neff_;
    neff_i_ = rhs.neff_i_;
    neff_d_ = rhs.neff_d_;
  }

  Real GetProb(HMMTransition transition) const {
    return trans_[transition];
  }

  void SetProb(HMMTransition transition, Real prob) {
    trans_[transition] = prob;
  }

  Real GetNeff() const {
    return neff_;
  }

  void SetNeff(Real neff) {
    neff_ = neff;
  }

  Real GetNeff_I() const {
    return neff_i_;
  }

  void SetNeff_I(Real neff) {
    neff_i_ = neff;
  }

  Real GetNeff_D() const {
    return neff_d_;
  }

  void SetNeff_D(Real neff) {
    neff_d_ = neff;
  }

  bool operator==(const HMMData& rhs) const {
    return (!memcmp(trans_, rhs.trans_, sizeof(trans_)) &&
            neff_ == rhs.neff_ &&
            neff_i_ == rhs.neff_i_ &&
            neff_d_ == rhs.neff_d_);
  }
  bool operator!=(const HMMData& rhs) const { return !(rhs == (*this)); }

  HMMData& operator=(const HMMData& rhs) {
    memcpy(trans_, rhs.trans_, 7*sizeof(Real));
    neff_ = rhs.neff_;
    neff_i_ = rhs.neff_i_;
    neff_d_ = rhs.neff_d_;
    return *this;
  }

  friend std::ofstream& operator<<(std::ofstream& os, HMMData& dat) {

    int16_t p_data[7];
    for (uint i = 0; i < 7; ++i) {
      p_data[i] = static_cast<int16_t>(dat.trans_[i]*10000);
    }
    os.write(reinterpret_cast<char*>(p_data), 7*sizeof(int16_t));

    float neff_data[3];
    neff_data[0] = dat.neff_;
    neff_data[1] = dat.neff_i_;
    neff_data[2] = dat.neff_d_;
    os.write(reinterpret_cast<char*>(neff_data), 3*sizeof(float));

    return os;
  }

  friend std::ifstream& operator>>(std::ifstream& is, HMMData& dat) {

    int16_t p_data[7];
    is.read(reinterpret_cast<char*>(p_data), 7*sizeof(int16_t));
    for (uint i = 0; i < 7; ++i) {
      dat.trans_[i] = p_data[i] * 0.0001;
    }

    float neff_data[3];
    is.read(reinterpret_cast<char*>(neff_data), 3*sizeof(float));
    dat.neff_ = neff_data[0];
    dat.neff_i_ = neff_data[1];
    dat.neff_d_ = neff_data[2];

    return is;
  }

private:
  Real trans_[7];
  Real neff_;
  Real neff_i_;
  Real neff_d_;
};

/// \brief Defines profile of 20 frequencies for one residue.
///
/// Frequencies are identified by the one-letter-code for that amino acid.
/// (possible codes: ACDEFGHIKLMNPQRSTVWY)
class DLLEXPORT_OST_SEQ ProfileColumn {
public:

  /// \brief Construct a profile with all frequencies set to 0.
  ProfileColumn() {
    memset(freq_, 0, sizeof(freq_));
  }

  ProfileColumn(const ProfileColumn& rhs) {
    memcpy(freq_, rhs.freq_, sizeof(freq_));
    if(rhs.hmm_data_) {
      // do deep copy
      hmm_data_ = HMMDataPtr(new HMMData(*rhs.hmm_data_));
    }
  }
  ProfileColumn& operator= (const ProfileColumn& rhs) {
    memcpy(freq_, rhs.freq_, sizeof(freq_));
    if(rhs.hmm_data_) {
      // do deep copy
      hmm_data_ = HMMDataPtr(new HMMData(*rhs.hmm_data_));
    } else if(hmm_data_) {
      hmm_data_ = HMMDataPtr();
    }
    return *this;
  }

  static ProfileColumn BLOSUMNullModel();
  static ProfileColumn HHblitsNullModel();

  void SetHMMData(HMMDataPtr p) {
    hmm_data_ = p;
  }

  HMMDataPtr GetHMMData() const{
    if(!hmm_data_) {
      throw Error("ProfileColumn has no HMM data set!");
    }
    return hmm_data_;
  }

  Real GetTransProb(HMMTransition transition) const {
    if(!hmm_data_) {
      throw Error("ProfileColumn has no HMM data set!");
    }
    return hmm_data_->GetProb(transition);    
  }

  /// \brief Translate one-letter-code to index (0-indexing).
  static int GetIndex(char ch);

  Real GetFreq(char ch) const;

  void SetFreq(char ch, Real freq);

  bool operator==(const ProfileColumn& rhs) const {
    return !memcmp(freq_, rhs.freq_, sizeof(freq_));
  }
  bool operator!=(const ProfileColumn& rhs) const { return !(rhs == (*this)); }

  Real* freqs_begin() { return freq_; }
  Real* freqs_end() { return freq_ + 20; }
  const Real* freqs_begin() const { return freq_; }
  const Real* freqs_end() const { return freq_ + 20; }

  /// \brief Get entropy for this column.
  Real GetEntropy() const;

  /// \brief Get column score as in Soeding-2005
  Real GetScore(const ProfileColumn& other,
                const ProfileColumn& null_model) const;

  // functions to feed streams with limited accuracy of internal data
  // not intended for python export

  friend std::ofstream& operator<<(std::ofstream& os, ProfileColumn& col) {
    int16_t data[20];
    //transform aa_freq
    for (uint i = 0; i < 20; ++i) {
      data[i] = static_cast<int16_t>(col.freq_[i]*10000);
    }
    os.write(reinterpret_cast<char*>(data), sizeof(data));
 
    if(col.hmm_data_) {
      bool has_hmm_data = true;
      os.write(reinterpret_cast<char*>(&has_hmm_data), 1);
      os << *col.hmm_data_;
    } else {
      bool has_hmm_data = false;
      os.write(reinterpret_cast<char*>(&has_hmm_data), 1);      
    }

    return os;
  }

  friend std::ifstream& operator>>(std::ifstream& is, ProfileColumn& col) {
    int16_t data[20];
    is.read(reinterpret_cast<char*>(data), sizeof(data));
    //transform aa_freq
    for (uint i = 0; i < 20; ++i) {
      col.freq_[i] = data[i] * 0.0001;
    }

    bool has_hmm_data;
    is.read(reinterpret_cast<char*>(&has_hmm_data), 1);
    if(has_hmm_data) {
      is >> *col.hmm_data_;
    }

    return is;
  }

private:
  Real freq_[20];
  HMMDataPtr hmm_data_;
};

/// \brief Provides a profile for a sequence.
///
/// Properties:
/// - consists of N ProfileColumn (N = #residues in sequence)
/// - contains sequence (String) of length N
/// - user must enforce consistency between sequence length and #columns
/// - contains a null_model to use for this sequence
class DLLEXPORT_OST_SEQ ProfileHandle { 
public:
  /// \brief Constructs an empty profile handle (sequence = '', 0 columns).
  ProfileHandle(): null_model_(ProfileColumn::HHblitsNullModel()), neff_(1.0) {}

  // uses compiler-generated copy- and assignment operators (work here!)

  const std::vector<ProfileColumn>& GetColumns() const { return columns_; }

  const ProfileColumn& GetNullModel() const { return null_model_; } 

  void SetNullModel(const ProfileColumn& null_model) {
    null_model_ = null_model;
  }

  String GetSequence() const { return seq_; }

  void SetSequence(const String& seq) {
    if (seq.length() != columns_.size()) {
      throw Error("ProfileHandle - Inconsistency between number of columns and "
                  " seq. length.");
    }
    seq_ = seq;
  }

  Real GetNeff() const { return neff_; }

  void SetNeff(Real neff) { neff_ = neff; }

  /// \brief Extract subset of profile for columns from until to-1 (0-indexing).
  /// Null model is copied from this profile.
  /// \throw Error if to <= from or to > size().
  ProfileHandlePtr Extract(uint from, uint to) const;

  /// \brief Compute average entropy over all columns.
  Real GetAverageEntropy() const;

  /// \brief Compute score comparing columns other[i] and this->at(i+offset)
  /// Column score as in Soeding-2005, null model of this object used, 
  /// result normalized by other.size()
  Real GetAverageScore(const ProfileHandle& other, uint offset = 0) const;

  // \brief Can only add column with an associated olc
  void AddColumn(const ProfileColumn& c, char olc='X') {
    columns_.push_back(c);
    seq_ += olc;
  } 

  // some functions to make it behave like a vector

  void clear() { seq_ = ""; columns_.clear(); }

  size_t size() const { return columns_.size(); }

  bool empty() const { return columns_.empty(); }

  ProfileColumn& operator[](size_t index) { return columns_[index]; }

  const ProfileColumn& operator[](size_t index) const { return columns_[index]; }

  ProfileColumn& at(size_t index) { return columns_.at(index); }

  const ProfileColumn& at(size_t index) const { return columns_.at(index); }

  ProfileColumn& back() { return columns_.back(); }

  const ProfileColumn& back() const { return columns_.back(); }

  bool operator==(const ProfileHandle& other) const {
    return seq_ == other.seq_ &&
           columns_ == other.columns_ &&
           null_model_ == other.null_model_;
  }
  
  bool operator!=(const ProfileHandle& other) const {
    return !(other == (*this));
  }

  ProfileColumnList::iterator columns_begin() { return columns_.begin(); }
  ProfileColumnList::iterator columns_end() { return columns_.end(); }
  ProfileColumnList::const_iterator columns_begin() const {
    return columns_.begin();
  }
  ProfileColumnList::const_iterator columns_end() const {
    return columns_.end();
  }

  // functions to feed streams with limited accuracy of internal data
  // not intended for python export

  friend std::ofstream& operator<<(std::ofstream& os, ProfileHandle& prof) {
    // null model
    os << prof.null_model_;
    // num. columns/residues
    uint32_t size = prof.size();
    os.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    // sequence
    if (prof.seq_.length() != size) {
      throw Error("ProfileHandle - Inconsistency between number of columns and "
                  " seq. length.");
    }
    os.write(prof.seq_.c_str(), size);
    // columns
    for(uint i = 0; i < size; ++i){
      os << prof.columns_[i];
    }

    return os;
  }

  friend std::ifstream& operator>>(std::ifstream& is, ProfileHandle& prof) {
    // null model
    is >> prof.null_model_;
    // num. columns/residues
    uint32_t size;
    is.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    // sequence
    std::vector<char> tmp_buf(size);
    is.read(&tmp_buf[0], size);
    prof.seq_.assign(&tmp_buf[0], size);
    // columns
    prof.columns_.resize(size);
    for(uint i = 0; i < size; ++i){
      is >> prof.columns_[i];
    }

    return is;
  }

private:
  String             seq_;
  ProfileColumn      null_model_; 
  ProfileColumnList  columns_;
  Real               neff_;
};

/// \brief Contains a DB of profiles (identified by a unique name (String)).
class DLLEXPORT_OST_SEQ ProfileDB {
public:
  /// \brief Saves all profiles in DB with limited accuracy of internal data.
  /// Binary format with fixed-width integers (should be portable).
  void Save(const String& filename) const;

  static ProfileDBPtr Load(const String& filename);

  void AddProfile(const String& name, ProfileHandlePtr prof);

  ProfileHandlePtr GetProfile(const String& name) const;

  size_t size() const { return data_.size(); }

  std::vector<String> GetNames() const;

private:
  std::map<String, ProfileHandlePtr> data_;
};

}}

#endif
