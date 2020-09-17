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

#ifndef OST_SEQ_ALG_HMM_PSEUDO_COUNTS_HH
#define OST_SEQ_ALG_HMM_PSEUDO_COUNTS_HH

#include <ost/seq/profile_handle.hh>

namespace ost{ namespace seq{ namespace alg{


class ContextProfileDB;
typedef boost::shared_ptr<ContextProfileDB> ContextProfileDBPtr;


class ContextProfile{

public:

ContextProfile(int length): length_(length), 
                            data_(ContextProfile::DataSize(length), 0.0) { }

ContextProfile(int length, Real* data): length_(length), 
                            data_(ContextProfile::DataSize(length), 0.0) {
  memcpy(&data_[0], data, data_.size() * sizeof(Real));
}

void SetWeight(int pos, char olc, Real weight) {
  if(pos >= length_) {
  	throw Error("Tried to access invalid pos in ContextProfile");
  }
  int olc_idx = ProfileColumn::GetIndex(olc);
  if(olc_idx != -1) {
    data_[pos*20 + olc_idx] = weight;
  } else {
  	throw Error("Invalid one letter code in ContextProfile");
  }
}

void SetPseudoCount(char olc, Real count) {
  int olc_idx = ProfileColumn::GetIndex(olc);
  if(olc_idx != -1) {
    data_[length_*20 + olc_idx] = count;
  } else {
  	throw Error("Invalid one letter code in ContextProfile");
  }
}

void SetBias(Real bias) { data_.back() = bias; }

const Real* GetWeights(int pos) const{ 
  if(pos >= length_) {
  	throw Error("Tried to access invalid pos in ContextProfile");
  }
  return &data_[pos*20];
}

Real GetWeight(int pos, char olc) {
  if(pos >= length_) {
    throw Error("Tried to access invalid pos in ContextProfile");
  }
  int olc_idx = ProfileColumn::GetIndex(olc);
  if(olc_idx != -1) {
    return data_[pos*20 + olc_idx];
  } else {
    throw Error("Invalid one letter code in ContextProfile");
  }  
}

const Real* GetPseudoCounts() const { return &data_[length_*20]; }

Real GetPseudoCount(char olc) {
  int olc_idx = ProfileColumn::GetIndex(olc);
  if(olc_idx != -1) {
    return data_[length_*20 + olc_idx];
  } else {
    throw Error("Invalid one letter code in ContextProfile");
  }    
}

Real GetBias() const { return data_.back(); }

const std::vector<Real>& GetData() const { return data_; }

int GetLength() const { return length_; }

static int DataSize(int length) { return (length+1)*20+1; }

private:
int length_;
// data organisation:
// context weights in chunks of 20 (length_ chunks)
// followed by 20 elements representing the context pseudo counts
// last element is the bias
std::vector<Real> data_;
};


class ContextProfileDB {

public:

ContextProfileDB() { } 

void Save(const String& filename) const;

static ContextProfileDBPtr Load(const String& filename);

static ContextProfileDBPtr FromCRF(const String& filename);

void AddProfile(const ContextProfile& profile){

  // enforce same length for all profiles
  if(!profiles_.empty()) {
    if(profile.GetLength() != profiles_[0].GetLength()) {
      throw Error("Require all profiles to be of same length");
    }
  }
  profiles_.push_back(profile);
}

const ContextProfile& operator [](int idx) const {
  return profiles_[idx];
}

const ContextProfile& at(int idx) const {
  return profiles_.at(idx);
}

size_t size() const {
  return profiles_.size();
}

size_t profile_length() const {
  if(profiles_.empty()) {
    throw Error("DB must contain profiles to get profile length");
  }
  return profiles_[0].GetLength();
}

private:
std::vector<ContextProfile> profiles_;
};

void AddTransitionPseudoCounts(ost::seq::ProfileHandle& profile,
                               Real gapb = 1.0, Real gabd = 0.15, 
                               Real gape = 1.0);

void AddAAPseudoCounts(ost::seq::ProfileHandle& profile,
                       Real a = 1.0, Real b = 1.5, Real c = 1.0);

void AddAAPseudoCounts(ost::seq::ProfileHandle& profile, 
                       const ContextProfileDB& db,
                       Real a = 0.9, Real b = 4.0, Real c = 1.0);

void AddNullPseudoCounts(ost::seq::ProfileHandle& profile);

}}} // ns

#endif
