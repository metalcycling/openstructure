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
  Author: Marco Biasini, Juergen Haas, Gerardo Tauriello
 */

#ifndef OST_SEQ_ALG_DISTANCE_MAP_HH
#define OST_SEQ_ALG_DISTANCE_MAP_HH

#include <numeric>
#include <cmath>
#include <algorithm>

#include <ost/tri_matrix.hh>
#include <ost/integrity_error.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/module_config.hh>

namespace ost { namespace seq { namespace alg {

///////////////////////////////////////////////////////////////////////////////
// HELPERS
namespace {

template <typename Pair>
struct sum_first {
  typename Pair::first_type operator()(const typename Pair::first_type& x,
                                       const Pair& y) const
  {
    return x + y.first;
  }
};

template <typename Pair>
struct CompareFirst {
  bool operator()(const Pair& x, const Pair& y) const {
    return (x.first < y.first);
  }
};

} // anon ns
///////////////////////////////////////////////////////////////////////////////

/// \brief Container for a pair wise distance for each structure.
/// Each structure is identified by its index in the originally used alignment.
struct DLLEXPORT_OST_SEQ_ALG Distances {

  typedef std::vector<std::pair<Real,int> > FloatArray;

  uint GetDataSize() const { return values_.size(); }

  void Add(Real dist, int index) {
    values_.push_back(std::make_pair(dist,index));
  }

  Real GetAverage() const {
    if (values_.empty()) throw IntegrityError("List of distances empty!");
    Real sum = std::accumulate(values_.begin(), values_.end(), 0.0,
                               sum_first<FloatArray::value_type>());
    return sum/values_.size();
  }

  std::pair<Real,int> GetMin() const {
    if (values_.empty()) throw IntegrityError("List of distances empty!");
    return *std::min_element(values_.begin(), values_.end(),
                             CompareFirst<FloatArray::value_type>());
  }

  std::pair<Real,int> GetMax() const {
    if (values_.empty()) throw IntegrityError("List of distances empty!");
    return *std::max_element(values_.begin(), values_.end(),
                             CompareFirst<FloatArray::value_type>());
  }

  std::pair<Real,int> GetDataElement(uint index) const {
    if (values_.size() != 0) {
      if (values_.size() > index) {
        return values_[index];
      } else {
        std::stringstream error_message;
        error_message << "size of distances data vector: " << values_.size()
                      << " is smaller than requested: " << index+1;
        throw IntegrityError(error_message.str());
      }
    } else {
        throw IntegrityError("List of distances empty!");
    }
  }

  Real GetStdDev() const {
    if (values_.empty()) throw IntegrityError("List of distances empty!");
    const Real avg = this->GetAverage();
    return this->do_get_stddev(avg);
  }

  Real GetWeightedStdDev(Real sigma) const {
    if (values_.empty()) throw IntegrityError("List of distances empty!");
    const Real avg = this->GetAverage();
    const Real weight = std::exp(avg/(-2*sigma));
    return this->do_get_stddev(avg) * weight;
  }

  Real GetNormStdDev() const {
    if (values_.empty()) throw IntegrityError("List of distances empty!");
    const Real avg = this->GetAverage();
    return this->do_get_stddev(avg) / avg;
  }

  bool operator==(const Distances& rhs) const { return values_ == rhs.values_; }
  bool operator!=(const Distances& rhs) const { return values_ != rhs.values_; }

private:

  Real do_get_stddev(Real avg) const {
    Real var = 0.0;
    for (FloatArray::const_iterator i = values_.begin(),
         e = values_.end(); i != e; ++i) {
      const Real diff = i->first - avg;
      var += diff*diff;
    }
    return std::sqrt(var/values_.size());
  }

  FloatArray values_;
};

/// \brief Container for pair wise distances of a set of structures.
/// Essentially a symmetric GetSize x GetSize matrix containing up to
/// GetNumStructures distances (see Distances).
class DistanceMap;
typedef boost::shared_ptr<DistanceMap> DistanceMapPtr;

class DLLEXPORT_OST_SEQ_ALG DistanceMap {
public:
  DistanceMap(uint nresidues, uint num_structures)
              : dist_(nresidues), num_structures_(num_structures) { }

  void AddDistance(uint i, uint j, Real distance, int index) {
    dist_(i ,j).Add(distance,index);
  }

  const Distances& GetDistances(uint i, uint j) const {
    return dist_(i, j);
  }

  uint GetSize() const { return dist_.GetSize(); }

  uint GetNumStructures() const { return num_structures_; }

private:
  TriMatrix<Distances> dist_;
  uint num_structures_;
};


/// \brief create distance map from a multiple sequence alignment.
///
/// The algorithm requires that the sequence alignment consists of at least two
/// sequences. The sequence at index 0 serves as a frame of reference. All the
/// other sequences must have an attached view and a properly set
/// \ref seq::SequenceHandle::SetSequenceOffset() "sequence offset".
///
/// For each of the attached views, the C-alpha distance pairs are extracted and
/// mapped onto the corresponding C-alpha distances in the reference sequence.
DistanceMapPtr DLLEXPORT_OST_SEQ_ALG
CreateDistanceMap(const seq::AlignmentHandle& alignment);

}}}

#endif
