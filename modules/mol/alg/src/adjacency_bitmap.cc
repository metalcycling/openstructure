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
#include <ost/img/image_factory.hh>
#include <ost/mol/atom_view.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/residue_view.hh>
#include "adjacency_bitmap.hh"

namespace ost { namespace mol { namespace alg {

const uint8_t AdjacencyBitmap::NUMBER_OF_BITS_SET[] = {
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 2, 
  1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 
  1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 
  1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 
  1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 
  2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 3, 3, 3, 3, 4 
};

const uint8_t AdjacencyBitmap::NUMBER_OF_EDGES_SET[] = {
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 
  1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 
  1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 
  1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 
  1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 
  1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 
  1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 
  2, 2, 3, 3, 2, 2, 3, 3, 3, 3, 4, 4, 3, 3, 4, 4, 
  2, 2, 3, 3, 2, 2, 3, 3, 3, 3, 4, 4, 3, 3, 4, 4, 
  1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 
  1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 
  2, 2, 3, 3, 2, 2, 3, 3, 3, 3, 4, 4, 3, 3, 4, 4, 
  2, 2, 3, 3, 2, 2, 3, 3, 3, 3, 4, 4, 3, 3, 4, 4
};

OverlapResult AdjacencyBitmap::Overlap(size_t vert_i, size_t vert_j) const {
  // To determine the union of the direct neighbors of vert_i and vert_j,
  // we calculated the bitwise or of the edge_bits and combine it with
  // the union of the defined bits. For that purpose, we will need a mask
  // for the edge bits. In binary the mask looks like 0101010101.
  uint64_t edge_mask =  0x5555555555555555;
  uint64_t* pi = data_+storage_size_*vert_i; 
  uint64_t* ei = pi + storage_size_;
  uint64_t* pj = data_+storage_size_*vert_j; 
  size_t nominator = 0;
  size_t denominator = 0;
  size_t edges_set = 0;
  for (; pi != ei; ++pi, ++pj) {
    uint64_t intersection_ij = *pi &  *pj;
    uint64_t union_ij = ((*pi | *pj) & edge_mask) | intersection_ij;
    for (size_t i = 0; i < sizeof(uint64_t); ++i) {
      nominator += NUMBER_OF_BITS_SET[(intersection_ij >> 8*i) & 0xff];
      edges_set += NUMBER_OF_EDGES_SET[(intersection_ij >> 8*i) & 0xff];
      denominator += NUMBER_OF_BITS_SET[(union_ij >> 8*i) & 0xff];
    }
  }
  return OverlapResult(nominator, denominator, edges_set);
}

void AdjacencyBitmap::dump() const {

  for (size_t i=0; i<this->size(); i+=10) {
    for (size_t j=0;j<this->size(); j+=10) {
      if (this->defined(i, j)) {
        std::cout << (this->get(i, j) ? '*' : '.');
      } else {
        std::cout << 'u';
      }
    }
    std::cout << std::endl;
  }
}

}}}

