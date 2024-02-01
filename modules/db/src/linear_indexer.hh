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


#ifndef OST_DB_LINEAR_INDEXER_HH
#define OST_DB_LINEAR_INDEXER_HH

#include <fstream>
#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <ost/base.hh>
#include <ost/stdint.hh>


namespace ost { namespace db {

class LinearIndexer;
typedef boost::shared_ptr<LinearIndexer> LinearIndexerPtr;

struct AssemblyInfo {

  void Read(std::ifstream& in_stream);

  void Write(std::ofstream& out_stream) const;

  bool operator==(const AssemblyInfo& other) const {
    return (data_start == other.data_start &&
            chain_names == other.chain_names &&
            chain_lengths == other.chain_lengths);
  }

  bool operator!=(const AssemblyInfo& other) const {
    return !(*this==other);
  }

  uint64_t data_start;
  std::vector<String> chain_names;
  std::vector<uint> chain_lengths;
};

class LinearIndexer {

public:

  LinearIndexer(): current_num_residues_(0) { }

  static LinearIndexerPtr Load(const String& filename);

  void Save(const String& filename) const;

  void AddAssembly(const String& name,
                   const std::vector<String>& chain_names,
                   const std::vector<uint>& chain_lengths);

  void RemoveAssembly(const String& name);

  std::vector<String> GetAssemblies() const;

  std::vector<String> GetChainNames(const String& name) const;

  std::vector<uint> GetChainLengths(const String& name) const;

  std::pair<uint64_t, uint64_t> GetDataRange(const String& name) const;

  std::pair<uint64_t, uint64_t> GetDataRange(const String& name,
                                             const String& chain_name) const;

  uint64_t GetNumResidues() const { return current_num_residues_; }

  bool operator==(const LinearIndexer& other) const {
    return (assemblies_ == other.assemblies_ &&
            idx_mapper_ == other.idx_mapper_ &&
            current_num_residues_ == other.current_num_residues_);

  }

  bool operator!=(const LinearIndexer& other) const {
    return !(*this == other);
  }

private:

  std::vector<AssemblyInfo> assemblies_;
  std::map<String, int> idx_mapper_;
  uint64_t current_num_residues_;
};


}} // ns

#endif
