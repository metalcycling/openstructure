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


#include <ost/db/linear_indexer.hh>
#include <ost/message.hh>
#include <numeric>

namespace ost { namespace db {

void AssemblyInfo::Read(std::ifstream& in_stream) {
  uint16_t num;
  in_stream.read(reinterpret_cast<char*>(&num), sizeof(uint16_t));
  chain_names.resize(num);
  chain_lengths.resize(num);
  in_stream.read(reinterpret_cast<char*>(&data_start), sizeof(uint64_t));

  for(uint i = 0; i < num; ++i) {
    uint16_t size;
    in_stream.read(reinterpret_cast<char*>(&size), sizeof(uint16_t));
    chain_names[i].resize(size);
    in_stream.read(&chain_names[i][0], size);
    uint32_t chain_length;
    in_stream.read(reinterpret_cast<char*>(&chain_length), sizeof(uint32_t));
    chain_lengths[i] = chain_length;
  }
}


void AssemblyInfo::Write(std::ofstream& out_stream) const {
  uint16_t num = chain_names.size();
  out_stream.write(reinterpret_cast<const char*>(&num), sizeof(uint16_t));
  out_stream.write(reinterpret_cast<const char*>(&data_start), sizeof(uint64_t));

  for(uint i = 0; i < num; ++i) {
    uint16_t size = chain_names[i].size();
    out_stream.write(reinterpret_cast<const char*>(&size), sizeof(uint16_t));
    out_stream.write(&chain_names[i][0], size);
    uint32_t chain_length = chain_lengths[i];
    out_stream.write(reinterpret_cast<const char*>(&chain_length), sizeof(uint32_t));
  }
}

LinearIndexerPtr LinearIndexer::Load(const String& filename) {

  std::ifstream in_stream(filename.c_str(), std::ios::binary);
  if (!in_stream) {
    throw ost::Error("Could not open " + filename);
  }

  LinearIndexerPtr loaded_db(new LinearIndexer);

  in_stream.read(reinterpret_cast<char*>(&loaded_db->current_num_residues_), 
                 sizeof(uint64_t));

  uint32_t num_assemblies;
  in_stream.read(reinterpret_cast<char*>(&num_assemblies), sizeof(uint32_t));

  loaded_db->assemblies_.resize(num_assemblies);
  for(uint i = 0; i < num_assemblies; ++i) {
    loaded_db->assemblies_[i].Read(in_stream);
  }

  for(uint i = 0; i < num_assemblies; ++i) {
    uint16_t size;
    in_stream.read(reinterpret_cast<char*>(&size), sizeof(uint16_t));
    String assembly_name(size, 'X');
    in_stream.read(&assembly_name[0], size);
    uint32_t idx;
    in_stream.read(reinterpret_cast<char*>(&idx), sizeof(uint32_t));
    loaded_db->idx_mapper_[assembly_name] = idx;
  }

  return loaded_db;
}


void LinearIndexer::Save(const String& filename) const {

  std::ofstream out_stream(filename.c_str(), std::ios::binary);
  if (!out_stream) {
    throw ost::Error("Could not open " + filename);
  }

  out_stream.write(reinterpret_cast<const char*>(&current_num_residues_), 
                   sizeof(uint64_t));

  uint32_t num_assemblies = assemblies_.size();
  out_stream.write(reinterpret_cast<const char*>(&num_assemblies), 
                   sizeof(uint32_t));

  for(uint i = 0; i < num_assemblies; ++i) {
    assemblies_[i].Write(out_stream);
  }

  for(std::map<String, int>::const_iterator it = idx_mapper_.begin();
      it != idx_mapper_.end(); ++it) {
    uint16_t size = it->first.size();
    out_stream.write(reinterpret_cast<const char*>(&size), sizeof(uint16_t));
    out_stream.write(&(it->first[0]), size);
    uint32_t idx = it->second;
    out_stream.write(reinterpret_cast<const char*>(&idx), sizeof(uint32_t));
  }
}


void LinearIndexer::AddAssembly(const String& name,
                                const std::vector<String>& chain_names,
                                const std::vector<uint>& chain_lengths) {

  if(chain_names.size() != chain_lengths.size()) {
    throw ost::Error("Chain names and lengths must be consistent in size!");
  }

  std::map<String, int>::iterator it = idx_mapper_.find(name);

  if(it != idx_mapper_.end()) {
    throw ost::Error("The assembly with name \"" + name +
                     "\" is already present in the indexer!");
  }
    
  // add new assembly info
  idx_mapper_[name] = assemblies_.size();
  assemblies_.push_back(AssemblyInfo());

  AssemblyInfo& assembly = assemblies_.back();
  assembly.data_start = current_num_residues_;
  assembly.chain_names = chain_names;
  assembly.chain_lengths = chain_lengths;
  current_num_residues_ += std::accumulate(chain_lengths.begin(), 
                                           chain_lengths.end(), 0);
}


void LinearIndexer::RemoveAssembly(const String& name) {

  std::map<String, int>::iterator idx_it = idx_mapper_.find(name);

  if(idx_it == idx_mapper_.end()) {
    throw ost::Error("The assembly with name \"" + name +
                     "\" is not present in the indexer!");
  }

  int assembly_idx = idx_it->second;
  uint64_t num_removed_res = 
  std::accumulate(assemblies_[assembly_idx].chain_lengths.begin(),
                  assemblies_[assembly_idx].chain_lengths.end(), 0);

  // Let's remove that assembly
  assemblies_.erase(assemblies_.begin() + assembly_idx);
  idx_mapper_.erase(idx_it);

  // Reduce all data start indices coming after the removed assembly
  for(uint i = assembly_idx; i < assemblies_.size(); ++i) {
    assemblies_[i].data_start -= num_removed_res;
  }

  // update the idx_mapper and numer of total residues
  for(std::map<String, int>::iterator it = idx_mapper_.begin();
      it != idx_mapper_.end(); ++it) {
    if(it->second > assembly_idx) {
      it->second -= 1;
    }
  }

  current_num_residues_ -= num_removed_res;
}


std::vector<String> LinearIndexer::GetAssemblies() const {

  std::vector<String> return_vec;
  for(std::map<String, int>::const_iterator it = idx_mapper_.begin();
      it != idx_mapper_.end(); ++it) {
    return_vec.push_back(it->first);
  }

  return return_vec;
}


std::vector<String> LinearIndexer::GetChainNames(const String& name) const {

  std::map<String, int>::const_iterator idx_it = idx_mapper_.find(name);

  if(idx_it == idx_mapper_.end()) {
    throw ost::Error("The assembly with name \"" + name +
                     "\" is not present in the indexer!");
  }

  const AssemblyInfo& assembly = assemblies_[idx_it->second];

  return assembly.chain_names;
}


std::vector<uint> LinearIndexer::GetChainLengths(const String& name) const {

  std::map<String, int>::const_iterator idx_it = idx_mapper_.find(name);

  if(idx_it == idx_mapper_.end()) {
    throw ost::Error("The assembly with name \"" + name +
                     "\" is not present in the indexer!");
  }

  const AssemblyInfo& assembly = assemblies_[idx_it->second];

  return assembly.chain_lengths;
}


std::pair<uint64_t, uint64_t> LinearIndexer::GetDataRange(
                                              const String& name) const {

  std::map<String, int>::const_iterator idx_it = idx_mapper_.find(name);

  if(idx_it == idx_mapper_.end()) {
    throw ost::Error("The assembly with name \"" + name +
                     "\" is not present in the indexer!");
  }

  const AssemblyInfo& assembly = assemblies_[idx_it->second];

  uint64_t num_res = std::accumulate(assembly.chain_lengths.begin(), 
                                     assembly.chain_lengths.end(), 0);

  return std::make_pair(assembly.data_start, assembly.data_start + num_res);
}


std::pair<uint64_t, uint64_t> LinearIndexer::GetDataRange(
                                              const String& name,
                                              const String& chain_name) const {

  std::map<String, int>::const_iterator idx_it = idx_mapper_.find(name);

  if(idx_it == idx_mapper_.end()) {
    throw ost::Error("The assembly with name \"" + name +
                     "\" is not present in the indexer!");
  }

  const AssemblyInfo& assembly = assemblies_[idx_it->second];

  int ch_idx = -1;
  int summed_length = 0;
  for(uint i = 0; i < assembly.chain_names.size(); ++i) {
    if(assembly.chain_names[i] == chain_name) { 
      ch_idx = i;
      break;
    }
    summed_length += assembly.chain_lengths[i];
  }

  if(ch_idx == -1) {
    throw ost::Error("Could not find chain with name \"" + chain_name +
                     "\" in specified assembly!");
  }

  uint64_t start = assembly.data_start + summed_length;
  uint64_t end = start + assembly.chain_lengths[ch_idx];

  return std::make_pair(start, end);    
}

}} //ns
