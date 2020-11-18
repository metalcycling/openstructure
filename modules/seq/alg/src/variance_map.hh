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

#ifndef OST_SEQ_ALG_VARIANCE_MAP_HH
#define OST_SEQ_ALG_VARIANCE_MAP_HH

#include <algorithm>

#include <ost/seq/alignment_handle.hh>
#include <ost/tri_matrix.hh>
#include <ost/integrity_error.hh>
#include <ost/seq/alg/module_config.hh>
#include <ost/seq/alg/distance_map.hh>

namespace ost { namespace seq { namespace alg {

class VarianceMap;
class Dist2Mean;
class MeanlDDT;
typedef boost::shared_ptr<VarianceMap> VarianceMapPtr;
typedef boost::shared_ptr<Dist2Mean> Dist2MeanPtr;
typedef boost::shared_ptr<MeanlDDT> MeanlDDTPtr;

/// \brief Container for variances for each entry in a distance map.
/// Main functionality: Get/Set, Min, Max, ExportXXX
/// Get/Set and GetSize taken from TriMatrix.
class DLLEXPORT_OST_SEQ_ALG VarianceMap: public TriMatrix<Real> {
public:
  // all values initialized to 0 in constructor!
  VarianceMap(int nresidues): TriMatrix<Real>(nresidues, 0) { }

  Real Min() {
    if (this->GetSize() == 0) throw IntegrityError("Matrix empty");
    std::vector<Real>& data = this->Data();
    return *std::min_element(data.begin(), data.end());
  }

  Real Max() {
    if (this->GetSize() == 0) throw IntegrityError("Matrix empty");
    std::vector<Real>& data = this->Data();
    return *std::max_element(data.begin(), data.end());
  }

  void ExportDat(const String& file_name);
  void ExportCsv(const String& file_name);
  void ExportJson(const String& file_name);
  String GetJsonString();
};

/// \brief Container for distances to mean for N structures.
/// Main functionality: Get/Set, ExportXXX
class DLLEXPORT_OST_SEQ_ALG Dist2Mean {
public:
  // all values initialized to 0 in constructor!
  Dist2Mean(uint num_residues, uint num_structures)
            : num_residues_(num_residues), num_structures_(num_structures)
            , values_(num_residues * num_structures, 0) { }

  void Set(uint i_res, uint i_str, Real val) {
    values_[GetIndex(i_res, i_str)] = val;
  }
  
  Real Get(uint i_res, uint i_str) const {
    return values_[GetIndex(i_res, i_str)];
  }
  
  Real& operator()(uint i_res, uint i_str) {
    return values_[GetIndex(i_res, i_str)];
  }
  Real operator()(uint i_res, uint i_str) const {
    return values_[GetIndex(i_res, i_str)];
  }
  
  std::vector<Real>& Data() { return values_; }

  uint GetNumResidues() const { return num_residues_; }
  uint GetNumStructures() const { return num_structures_; }

  void Add(uint i_res, uint i_str, Real val) {
    values_[GetIndex(i_res, i_str)] += val;
  }

  void DivideBy(Real val) {
    for (uint i = 0; i < values_.size(); ++i) values_[i] /= val;
  }

  void ExportDat(const String& file_name);
  void ExportCsv(const String& file_name);
  void ExportJson(const String& file_name);
  String GetJsonString();
  
private:
  uint GetIndex(uint i_res, uint i_str) const {
    assert(i_res < num_residues_);
    assert(i_str < num_structures_);
    return (i_res * num_structures_ + i_str);
  }

  uint num_residues_;
  uint num_structures_;
  std::vector<Real> values_;
};


class DLLEXPORT_OST_SEQ_ALG MeanlDDT {
public:
  // all values initialized to 0 in constructor!
  MeanlDDT(uint num_residues, uint num_structures)
           : num_residues_(num_residues), num_structures_(num_structures)
            , values_(num_residues * num_structures, 0) { }

  void Set(uint i_res, uint i_str, Real val) {
    values_[GetIndex(i_res, i_str)] = val;
  }
  
  Real Get(uint i_res, uint i_str) const {
    return values_[GetIndex(i_res, i_str)];
  }
  
  Real& operator()(uint i_res, uint i_str) {
    return values_[GetIndex(i_res, i_str)];
  }
  Real operator()(uint i_res, uint i_str) const {
    return values_[GetIndex(i_res, i_str)];
  }
  
  std::vector<Real>& Data() { return values_; }

  uint GetNumResidues() const { return num_residues_; }
  uint GetNumStructures() const { return num_structures_; }

  void ExportDat(const String& file_name);
  void ExportCsv(const String& file_name);
  void ExportJson(const String& file_name);
  String GetJsonString();
  
private:
  uint GetIndex(uint i_res, uint i_str) const {
    assert(i_res < num_residues_);
    assert(i_str < num_structures_);
    return (i_res * num_structures_ + i_str);
  }

  uint num_residues_;
  uint num_structures_;
  std::vector<Real> values_;
};

/// \returns Variance measure for each entry in dmap.
/// \param dmap Distance map as created with CreateDistanceMap.
/// \param sigma Used for weighting of variance measure
///              (see Distances::GetWeightedStdDev)
VarianceMapPtr DLLEXPORT_OST_SEQ_ALG
CreateVarianceMap(const DistanceMapPtr dmap, Real sigma=25);

/// \returns Distances to mean for each structure in dmap.
///          Structures are in the same order as passed when creating dmap.
/// \param dmap Distance map as created with CreateDistanceMap.
Dist2MeanPtr DLLEXPORT_OST_SEQ_ALG
CreateDist2Mean(const DistanceMapPtr dmap);

MeanlDDTPtr DLLEXPORT_OST_SEQ_ALG
CreateMeanlDDTHA(const DistanceMapPtr dmap);

}}}

#endif
