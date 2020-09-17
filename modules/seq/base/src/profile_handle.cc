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

#include <fstream>
#include <iostream>
#include <ost/seq/profile_handle.hh>

namespace ost { namespace seq{ 

//////////////////////////////////////////////////////////////////////////////
// ProfileColumn
//////////////////////////////////////////////////////////////////////////////

ProfileColumn ProfileColumn::BLOSUMNullModel() {
  ProfileColumn col;

  col.freq_[0] = 0.0732;
  col.freq_[1] = 0.0250;
  col.freq_[2] = 0.0542;
  col.freq_[3] = 0.0547;
  col.freq_[4] = 0.0447;
  col.freq_[5] = 0.0751;
  col.freq_[6] = 0.0261;
  col.freq_[7] = 0.0700;
  col.freq_[8] = 0.1011;
  col.freq_[9] = 0.0584;
  col.freq_[10] = 0.0246;
  col.freq_[11] = 0.0463;
  col.freq_[12] = 0.0351;
  col.freq_[13] = 0.0326;
  col.freq_[14] = 0.0484;
  col.freq_[15] = 0.0573;
  col.freq_[16] = 0.0505;
  col.freq_[17] = 0.0758;
  col.freq_[18] = 0.0124;
  col.freq_[19] = 0.0345;

  return col;
}

ProfileColumn ProfileColumn::HHblitsNullModel() {
  ProfileColumn col;

  col.freq_[0] = 0.07662717998027802;
  col.freq_[1] = 0.01886688359081745;
  col.freq_[2] = 0.05399613827466965;
  col.freq_[3] = 0.05978801101446152;
  col.freq_[4] = 0.03493943437933922;
  col.freq_[5] = 0.0754152461886406;
  col.freq_[6] = 0.03682935610413551;
  col.freq_[7] = 0.050485048443078995;
  col.freq_[8] = 0.05958116054534912;
  col.freq_[9] = 0.09992572665214539;
  col.freq_[10] = 0.021959668025374413;
  col.freq_[11] = 0.040107060223817825;
  col.freq_[12] = 0.045310840010643005;
  col.freq_[13] = 0.03264486789703369;
  col.freq_[14] = 0.05129634961485863;
  col.freq_[15] = 0.04661700129508972;
  col.freq_[16] = 0.07105106115341187;
  col.freq_[17] = 0.07264462858438492;
  col.freq_[18] = 0.012473411858081818;
  col.freq_[19] = 0.039418045431375504;

  return col;
}

int ProfileColumn::GetIndex(char olc) {
  if (olc == 'A') return 0;
  if (olc >= 'C' && olc <= 'I') return (olc-'A') - 1;
  if (olc >= 'K' && olc <= 'N') return (olc-'A') - 2;
  if (olc >= 'P' && olc <= 'T') return (olc-'A') - 3;
  if (olc == 'V' ) return 17;
  if (olc == 'W' ) return 18;
  if (olc == 'Y' ) return 19;
  return -1;
}

Real ProfileColumn::GetFreq(char ch) const{
  int idx = ProfileColumn::GetIndex(ch);
  if (idx == -1) {
    throw Error("Invalid One Letter Code observed when getting frequency in "
                "ProfileColumn!");
  }
  return freq_[idx];
}

void ProfileColumn::SetFreq(char ch, Real freq){
  int idx = ProfileColumn::GetIndex(ch);
  if (idx == -1) {
    throw Error("Invalid One Letter Code observed when setting frequency in "
                "ProfileColumn!");
  }
  freq_[idx] = freq;
}

Real ProfileColumn::GetEntropy() const {
  Real entropy = 0.0;
  for (uint i = 0; i < 20; ++i) {
    if (freq_[i] > 0.0) {
      entropy -= log(freq_[i]) * freq_[i];
    }
  }
  return entropy;
}

Real ProfileColumn::GetScore(const ProfileColumn& other,
                             const ProfileColumn& null_model) const {
  Real summed_col_score = 0.001; // to avoid zero score
  for (uint i = 0; i < 20; ++i) {
    summed_col_score += freq_[i] * other.freq_[i] / null_model.freq_[i];
  }
  return std::log(summed_col_score);
}

//////////////////////////////////////////////////////////////////////////////
// ProfileHandle
//////////////////////////////////////////////////////////////////////////////

ProfileHandlePtr ProfileHandle::Extract(uint from, uint to) const {
  // check
  if (to <= from) {
    throw Error("Second index must be bigger than first one!");
  }
  if (to > this->size()) {
    throw Error("Invalid index!");
  }

  // get subset to return
  ProfileHandlePtr return_prof(new ProfileHandle);
  return_prof->SetNullModel(null_model_);
  return_prof->SetNeff(neff_);
  for (uint i = from; i < to; ++i) {
    return_prof->AddColumn(columns_[i], seq_[i]);
  }

  return return_prof;
}

Real ProfileHandle::GetAverageEntropy() const {
  Real n_eff = 0.0;
  int n = 0;
  for (std::vector<ProfileColumn>::const_iterator
       i = this->columns_begin(), e = this->columns_end(); i != e; ++i) {
    n += 1;
    n_eff += i->GetEntropy();
  }
  return (n > 0) ? n_eff/n : 0.0;
}

Real ProfileHandle::GetAverageScore(const ProfileHandle& other,
                                    uint offset) const {
  // check offsets
  const uint other_size = other.size();
  if (offset > this->size() || offset + other_size > this->size()) {
    throw Error("Invalid profile size / offset when computing score!");
  }
  // sum them up
  Real sum = 0;
  for (uint i = 0; i < other_size; ++i) {
    sum += columns_[offset + i].GetScore(other[i], null_model_);
  }
  if (other_size > 0) return sum / other_size;
  else                return sum;
}

//////////////////////////////////////////////////////////////////////////////
// ProfileDB
//////////////////////////////////////////////////////////////////////////////


void ProfileDB::Save(const String& filename) const{

  std::ofstream out_stream(filename.c_str(), std::ios::binary);

  uint32_t magic_number = 42;
  out_stream.write(reinterpret_cast<char*>(&magic_number), sizeof(uint32_t));

  uint8_t version = 1;
  out_stream.write(reinterpret_cast<char*>(&version), sizeof(uint8_t));

  //write out total size
  uint32_t total_size = data_.size();
  out_stream.write(reinterpret_cast<char*>(&total_size),4);

  //write out the data elements
  char string_size;
  for (std::map<String, ProfileHandlePtr>::const_iterator i = data_.begin();
       i != data_.end(); ++i){
    string_size = static_cast<char>(i->first.size());
    out_stream.write(reinterpret_cast<char*>(&string_size),1);
    out_stream.write(i->first.c_str(),string_size);
    out_stream << *i->second;
  }
  out_stream.close();
}

ProfileDBPtr ProfileDB::Load(const String& filename){

  std::ifstream in_stream(filename.c_str(), std::ios::binary);
  if (!in_stream){
    std::stringstream ss;
    ss << "the file '" << filename << "' does not exist.";
    throw Error(ss.str());
  }

  ProfileDBPtr db(new ProfileDB);

  uint32_t magic_number;
  in_stream.read(reinterpret_cast<char*>(&magic_number), sizeof(uint32_t));

  if(magic_number != 42) {
    std::stringstream ss;
    ss << "Could not read magic number in " << filename<<". Either the file ";
    ss << "is corrupt, does not contain a ProfileDB or is of an old version ";
    ss << "which is not supported anymore.";
    throw Error(ss.str());
  }

  uint8_t version;
  in_stream.read(reinterpret_cast<char*>(&version), sizeof(uint8_t));
  if(version != 1) {
    std::stringstream ss;
    ss << "ProfileDB in " << filename << " is of version " << version;
    ss << " but only version 1 can be read.";
    throw Error(ss.str());
  }

  //read in the total size
  uint32_t total_size;
  in_stream.read(reinterpret_cast<char*>(&total_size),4);

  //read in the single profiles
  char string_size;
  for(uint i = 0; i < total_size; ++i){
    in_stream.read(&string_size,1);
    String name(string_size,'X');
    in_stream.read(&name[0],string_size);
    ProfileHandlePtr prof(new ProfileHandle);
    in_stream >> *prof;
    db->AddProfile(name, prof);
  }

  return db;
}

void ProfileDB::AddProfile(const String& name, ProfileHandlePtr prof) {
  if (name.size() > 255) {
    throw Error("Name of Profile must be smaller than 256!");
  }
  if (name.empty()) {
    throw Error("Name must not be empty!");
  }
  data_[name] = prof;
}

ProfileHandlePtr ProfileDB::GetProfile(const String& name) const{
  std::map<String,ProfileHandlePtr>::const_iterator i = data_.find(name);
  if (i == data_.end()) {
    std::stringstream ss;
    ss << "Profile database does not contain an entry with name ";
    ss << name << "!";
    throw Error(ss.str());
  }
  return i->second;
}

std::vector<String> ProfileDB::GetNames() const{
  std::vector<String> return_vec;
  return_vec.reserve(this->size());
  for (std::map<String, ProfileHandlePtr>::const_iterator i = data_.begin(); 
       i != data_.end(); ++i){
    return_vec.push_back(i->first);
  }
  return return_vec;
}


}} //ns
