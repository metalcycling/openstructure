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
  Author: Gerardo Tauriello, Juergen Haas
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

#include <ost/mol/mol.hh>

#include <ost/seq/alg/distance_map.hh>
#include <ost/seq/alg/variance_map.hh>

namespace ost { namespace seq { namespace alg {

///////////////////////////////////////////////////////////////////////////////
// HELPERS
namespace {
// dump stuff using () operator
template <typename T>
void DumpCsv(const String& file_name, const T& data, uint num_rows,
             uint num_cols) {
  // open file
  std::ofstream out_file(file_name.c_str());
  if (!out_file) throw Error("The file '" + file_name + "' cannot be opened.");
  // dump each row
  for (uint row = 0; row < num_rows; ++row) {
    for (uint col = 0; col < num_cols; ++col) {
      out_file << data(row, col);
      if (col < num_cols-1) out_file << ";";
    }
    out_file << std::endl;
  }
  out_file.close();
}

template <typename T>
String GetJson(const T& data, uint num_rows, uint num_cols) {
  std::ostringstream out_stream;
  out_stream << "[";
  for (uint row = 0; row < num_rows; ++row) {
    out_stream << "[";
    for (uint col = 0; col < num_cols; ++col) {
      out_stream << data(row, col);
      if (col < num_cols-1) out_stream << ",";
    }
    out_stream << "]";
    if (row < num_rows-1) out_stream << ",";
  }
  out_stream << "]";
  return out_stream.str();
}


void FilllDDTValues(const std::vector<Real>& dist_diff,
                    const std::vector<int>& ref_distances,
                    const std::vector<std::pair<int, int> >& idx_mapping,
                    std::vector<Real>& local_lddt) {

  int len = dist_diff.size();
  std::vector<int> counts(len, 0);
  std::vector<int> total_counts(len, 0);

  for(auto it = ref_distances.begin(); it != ref_distances.end(); ++it) {
    Real diff = dist_diff[*it];
    int fulfilled = 0;
    fulfilled += static_cast<int>(diff < static_cast<Real>(4.0));
    fulfilled += static_cast<int>(diff < static_cast<Real>(2.0));
    fulfilled += static_cast<int>(diff < static_cast<Real>(1.0));
    fulfilled += static_cast<int>(diff < static_cast<Real>(0.5));

    const std::pair<int, int>& d_indices = idx_mapping[*it];
    counts[d_indices.first] += fulfilled;
    counts[d_indices.second] += fulfilled;
    total_counts[d_indices.first] += 4;
    total_counts[d_indices.second] += 4;
  }

  local_lddt.assign(len, std::numeric_limits<Real>::quiet_NaN());

  for(int i = 0; i < len; ++i) {
    if(counts[i] > 0) {
      local_lddt[i] = static_cast<Real>(counts[i]) / total_counts[i];
    }
  }

}


void FilllDDTValues(const std::vector<Real>& dist_one, 
                    const std::vector<Real>& dist_two,
                    const std::vector<int>& below_15_one, 
                    const std::vector<int>& below_15_two, 
                    const std::vector<std::pair<int, int> >& idx_mapping,
                    std::vector<Real>& local_lddt_one, 
                    std::vector<Real>& local_lddt_two) {

  // absolute dist differences are the same for both
  std::vector<Real> dist_diff(dist_one.size());
  for(uint i = 0; i < dist_one.size(); ++i) {
    dist_diff[i] = std::abs(dist_one[i] - dist_two[i]);
  }

  // estimate lDDT of structure 1 with reference distances (below 15) from 
  // structure 2 and vice versa
  FilllDDTValues(dist_diff, below_15_two, idx_mapping, local_lddt_one);
  FilllDDTValues(dist_diff, below_15_one, idx_mapping, local_lddt_two);
}

} // anon ns

///////////////////////////////////////////////////////////////////////////////
// IO
void VarianceMap::ExportDat(const String& file_name) {
  uint len = this->GetSize();
  if (len == 0) throw IntegrityError("Matrix empty");
  // dump it
  std::ofstream out_file(file_name.c_str());
  if (!out_file) throw Error("The file '" + file_name + "' cannot be opened.");
  for (uint rows = 0; rows < len; ++rows) {
    for (uint cols = 0; cols < len; ++cols) {
      out_file << (rows+1) << " " << (cols+1) << " " << this->Get(rows, cols)
               << std::endl;
    }
  }
  out_file.close();
}

void VarianceMap::ExportCsv(const String& file_name) {
  uint len = this->GetSize();
  if (len == 0) throw IntegrityError("Matrix empty");
  DumpCsv(file_name, *this, len, len);
}

void VarianceMap::ExportJson(const String& file_name) {
  if (this->GetSize() == 0) throw IntegrityError("Matrix empty");
  std::ofstream out_file(file_name.c_str());
  if (!out_file) throw Error("The file '" + file_name + "' cannot be opened.");
  out_file << this->GetJsonString() << std::endl;
  out_file.close();
}

String VarianceMap::GetJsonString() {
  uint len = this->GetSize();
  if (len == 0) throw IntegrityError("Matrix empty");
  return GetJson(*this, len, len);
}

void Dist2Mean::ExportDat(const String& file_name) {
  if (values_.size() == 0) throw IntegrityError("Matrix empty");
  // dump it
  std::ofstream out_file(file_name.c_str());
  if (!out_file) throw Error("The file '" + file_name + "' cannot be opened.");
  for (uint i_res = 0; i_res < num_residues_; ++i_res) {
    out_file << (i_res+1);
    for (uint i_str = 0; i_str < num_structures_; ++i_str) {
      out_file << " " << this->Get(i_res, i_str);
    }
    out_file << std::endl;
  }
  out_file.close();
}

void Dist2Mean::ExportCsv(const String& file_name) {
  if (values_.size() == 0) throw IntegrityError("Matrix empty");
  DumpCsv(file_name, *this, num_residues_, num_structures_);
}

void Dist2Mean::ExportJson(const String& file_name) {
  if (values_.size() == 0) throw IntegrityError("Matrix empty");
  std::ofstream out_file(file_name.c_str());
  if (!out_file) throw Error("The file '" + file_name + "' cannot be opened.");
  out_file << this->GetJsonString() << std::endl;
  out_file.close();
}

String Dist2Mean::GetJsonString() {
  if (values_.size() == 0) throw IntegrityError("Matrix empty");
  return GetJson(*this, num_residues_, num_structures_);
}

void MeanlDDT::ExportDat(const String& file_name) {
  if (values_.size() == 0) throw IntegrityError("Matrix empty");
  // dump it
  std::ofstream out_file(file_name.c_str());
  if (!out_file) throw Error("The file '" + file_name + "' cannot be opened.");
  for (uint i_res = 0; i_res < num_residues_; ++i_res) {
    out_file << (i_res+1);
    for (uint i_str = 0; i_str < num_structures_; ++i_str) {
      out_file << " " << this->Get(i_res, i_str);
    }
    out_file << std::endl;
  }
  out_file.close();
}

void MeanlDDT::ExportCsv(const String& file_name) {
  if (values_.size() == 0) throw IntegrityError("Matrix empty");
  DumpCsv(file_name, *this, num_residues_, num_structures_);
}

void MeanlDDT::ExportJson(const String& file_name) {
  if (values_.size() == 0) throw IntegrityError("Matrix empty");
  std::ofstream out_file(file_name.c_str());
  if (!out_file) throw Error("The file '" + file_name + "' cannot be opened.");
  out_file << this->GetJsonString() << std::endl;
  out_file.close();
}

String MeanlDDT::GetJsonString() {
  if (values_.size() == 0) throw IntegrityError("Matrix empty");
  return GetJson(*this, num_residues_, num_structures_);
}

///////////////////////////////////////////////////////////////////////////////
// Create maps
VarianceMapPtr CreateVarianceMap(const DistanceMapPtr dmap, Real sigma) {
  // setup
  uint len = dmap->GetSize();
  if (len == 0) {
    throw IntegrityError("empty distance map provided");
  }
  // get weighted std for each pair (symmetric storage!)
  VarianceMapPtr vmap(new VarianceMap(len));
  for (uint i_res = 0; i_res < len; ++i_res) {
    for (uint i_res2 = 0; i_res2 <= i_res; ++i_res2) {
      const Distances& di = dmap->GetDistances(i_res, i_res2);
      if (di.GetDataSize() > 0) {
        vmap->Set(i_res, i_res2, di.GetWeightedStdDev(sigma));
      } else {
        vmap->Set(i_res, i_res2, 0);
      }
    }
  }
  return vmap;
}

Dist2MeanPtr CreateDist2Mean(const DistanceMapPtr dmap) {
  // setup/check
  uint nstruc = dmap->GetNumStructures();
  uint len = dmap->GetSize();
  if (len == 0 || nstruc == 0) {
    throw IntegrityError("empty distance map provided");
  }
  // sum up all distances to mean for each residue (symmetric dmap!)
  Dist2MeanPtr dist2mean(new Dist2Mean(len, nstruc));
  for (uint i_res = 0; i_res < len; ++i_res) {
    for (uint i_res2 = 0; i_res2 <= i_res; ++i_res2) {
      const Distances& di = dmap->GetDistances(i_res, i_res2);
      const uint size = di.GetDataSize();
      if (size >= 2) {
        const Real avg = di.GetAverage();
        for (uint h = 0; h < size; ++h) {
          const std::pair<Real, int> ret = di.GetDataElement(h);
          const Real val = std::abs(ret.first - avg);
          dist2mean->Add(i_res, ret.second-1, val);
          if (i_res != i_res2) dist2mean->Add(i_res2, ret.second-1, val);
        }
      }
    }
  }
  // normalize
  dist2mean->DivideBy(len);
  return dist2mean;
}


MeanlDDTPtr CreateMeanlDDTHA(const DistanceMapPtr dmap) {
  // setup/check
  uint nstruc = dmap->GetNumStructures();
  uint len = dmap->GetSize();
  if (len <= 1 || nstruc == 0) {
    throw IntegrityError("empty distance map provided");
  }

  // the relevant pairwise distances are the off-diagonal elements from the
  // pairwise distance matrices -> the upper half without the diagonal elements
  uint n_off_diagonal = (len*len-len)/2;
  // store that information in a linear layout for each structure
  std::vector<std::vector<Real> > dist_data(dmap->GetNumStructures());
  // keep track of distances <= 15 for lDDT calculation. Indices stored in here
  // relate to the linear representations in dist_data
  std::vector<std::vector<int> > below_15(dmap->GetNumStructures());
  for(uint s_idx = 0; s_idx < nstruc; ++s_idx) {
    // the default value is simply a very large number this will trigger a 
    // distance difference above any relevant threshold if not set
    dist_data[s_idx].assign(n_off_diagonal, 10000000.0);
  }
  // keep track of which two residues are involved for each element in dist_data
  std::vector<std::pair<int,int> > off_diagonal_mapper(n_off_diagonal);

  uint off_diag_idx = 0;
  for (uint i_res = 0; i_res < len; ++i_res) {
    for (uint j_res = i_res + 1; j_res < len; ++j_res, ++off_diag_idx) {
      off_diagonal_mapper[off_diag_idx] = std::make_pair(i_res, j_res);
      const Distances& di = dmap->GetDistances(i_res, j_res);
      for (uint k = 0; k < di.GetDataSize(); ++k) {
        const std::pair<Real, int>& ret = di.GetDataElement(k);
        dist_data[ret.second-1][off_diag_idx] = ret.first;
        if(ret.first <= 15.0) {
          below_15[ret.second-1].push_back(off_diag_idx);
        }
      }
    }
  }

  std::vector<std::vector<Real> > values(nstruc);
  std::vector<std::vector<int> > counts(nstruc);
  for(uint i = 0; i < nstruc; ++i) {
    values[i].assign(len, 0);
    counts[i].assign(len, 0);
  }

  std::vector<Real> lddt_i, lddt_j;
  for(uint i = 0; i < nstruc; ++i) {
    for(uint j = i+1; j < nstruc; ++j) {
      FilllDDTValues(dist_data[i], dist_data[j], below_15[i], below_15[j], 
                     off_diagonal_mapper, lddt_i, lddt_j);
      for(uint k = 0; k < len; ++k) {
        if(!std::isnan(lddt_i[k])) {
          values[i][k] += lddt_i[k];
          counts[i][k] += 1;
        }
        if(!std::isnan(lddt_j[k])) {
          values[j][k] += lddt_j[k];
          counts[j][k] += 1;
        }
      }
    }
  }

  MeanlDDTPtr return_ptr(new MeanlDDT(len, nstruc));
  for(uint struc_idx = 0; struc_idx < nstruc; ++struc_idx) {
    for(uint res_idx = 0; res_idx < len; ++res_idx) {
      if(counts[struc_idx][res_idx] > 0) {
        return_ptr->Set(res_idx, struc_idx, 
                        values[struc_idx][res_idx]/counts[struc_idx][res_idx]);
      }
    }
  }

  return return_ptr;
}


 
}}}
